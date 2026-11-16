import json

import pandas as pd

from sources import normalize

DIAGNOSIS_OVERRIDES = {
    'CHUNG-JANSEN SYNDROME; CHUJANS': (['OMIM:617991', 'ORPHA:589905'], 'Curated PHIP synonym; NCBI GTR C4693860 verifies identity'),
    'FBXW7': (['OMIM:620012'], 'Gene-labelled cohort; mapped to FBXW7-related developmental delay/hypotonia/impaired language, not generic ORPHA intellectual disability'),
    'Baraitser-Winter syndrome': (['ORPHA:2995'], 'Broad cerebrofrontofacial syndrome; do not pick one ACTB/ACTG1 subtype'),
    'Hutchinson-Gilford progeria': (['OMIM:176670', 'ORPHA:740'], 'Exact established syndrome alias'),
    'Hyperphosphatasia with mental retardation syndrome': (['ORPHA:247262'], 'Mabry/hyperphosphatasia-intellectual disability syndrome umbrella; not a single molecular subtype'),
    'INTELLECTUAL DEVELOPMENTAL DISORDER, X-LINKED, SYNDROMIC, SNIJDERS BLOK TYPE; MRXSSB': (['OMIM:300958'], 'Exact Snijders Blok disease; spelling differs in annotation label'),
    'MANNOSIDOSIS, ALPHA B, LYSOSOMAL; MANSA': (['OMIM:248500', 'ORPHA:61'], 'Alpha-mannosidosis; explicitly reject beta-mannosidosis OMIM:248510'),
    'Joubert syndrome': (['ORPHA:140874'], 'Broad Joubert syndrome; no substitution of an individual numbered subtype'),
    'NEURODEVELOPMENTAL DISORDER WITH DYSMORPHIC FACIES AND DISTAL LIMB ANOMALIES; NEDDFL':
        (['OMIM:617755', 'ORPHA:686482'], 'Exact BPTF NEDDFL identity; current Orphadata alias adds BPTF-related prefix'),
    'Arthrogryposis, distal': (['ORPHA:97120'], 'Broad clinical group; subtype composition unknown'),
    'Cutis laxa': (['ORPHA:209'], 'Broad clinical group; subtype composition unknown'),
    'Mucopolysaccharidoses': (['ORPHA:79213'], 'Broad category; do not substitute MPS-II or another individual subtype'),
    'Autism, susceptiblity to': ([], 'Ambiguous susceptibility umbrella; no unique disease annotation or accepted typical facial phenotype assigned'),
}
UNRESOLVED_SCOPE = {'Arthrogryposis, distal', 'Cutis laxa', 'Mucopolysaccharidoses', 'Autism, susceptiblity to'}


def disease_mapping(syndromes, annotations, nomenclature):
    names = annotations[['disease_id', 'disease_name']].drop_duplicates()
    rows = []
    for disease in sorted(syndromes):
        display = disease.split('__')[0]
        stem = display.split(';')[0].strip()
        candidates = names.loc[names.disease_name.map(normalize).eq(normalize(stem)), 'disease_id'].tolist()
        note = 'Exact normalized clinical name only; punctuation/case/diacritic normalization, no fuzzy disease joins'
        if display in DIAGNOSIS_OVERRIDES:
            candidates, note = DIAGNOSIS_OVERRIDES[display]
        for record in nomenclature:
            if normalize(stem) in {normalize(t) for t in [record['name']] + record['synonyms']}:
                candidates = sorted(set(candidates + [record['orpha']]))
        record_names = names.loc[names.disease_id.isin(candidates)].drop_duplicates()
        annotation_count = int(annotations.disease_id.isin(candidates).sum())
        rows.append({'disease': disease, 'display_name': display, 'disease_ids': '|'.join(sorted(set(candidates))),
                     'annotation_disease_names': '|'.join(sorted(record_names.disease_name.unique())),
                     'annotation_rows': annotation_count, 'scope_unresolved': stem in UNRESOLVED_SCOPE,
                     'mapping_status': 'curated_identity' if candidates else 'unresolved_identity',
                     'scope_note': note, 'race_specific_clinical_source': False})
    return pd.DataFrame(rows)


def rules():
    rows = []

    def add(name, patterns, expected, related=(), note='', relation='normalized_anatomical_proxy'):
        rows.append({'family': name, 'patterns': patterns, 'expected': expected, 'related': related,
                     'note': note, 'relation': relation})

    add('eyebrow_polygon_area', ['^eb_thickness_(r|l|mean)$'], {'HP:0000574': 1, 'HP:0045074': -1, 'HP:0011229': 1},
        ['HP:0045075'], 'Polygon area is not hair density; sparse eyebrows are related but supply no geometry direction.')
    add('eyebrow_lateral_geometry', ['^eb_lateral_thickness_(r|l|mean)$'], {'HP:0007933': 1},
        ['HP:0005338', 'HP:0000574', 'HP:0045074'], 'Lateral landmark separation does not measure lateral hair abundance or whole-brow thickness.')
    add('eyebrow_arch', ['^eb_arch_(r|l|mean)$'], {'HP:0002553': 1, 'HP:0011228': -1},
        note='Arch height is divided by brow length; not a clinical diagnostic threshold.')
    add('eyebrow_spacing', ['^eb_inter_distance$'], {'HP:0000664': -1},
        note='Synophrys concerns medial eyebrow continuity; projected spacing is an incomplete proxy.')
    add('eyebrow_position', ['^eb_position_(r|l|mean)$'], {}, ['HP:0040296', 'HP:0010731'],
        'Brow-to-upper-lid distance changes with both eyebrow position and eyelid position; no unique expected direction.')
    add('eyebrow_medial_flare', ['^eb_medial_flare_(r|l|mean)$'], {}, ['HP:0010747'],
        'HPO describes medial hair configuration; landmark slope has no validated clinical direction mapping.')
    add('eyebrow_asymmetry', ['^eb_.*_asym$'], {'HP:6001499': 1},
        note='Ratio or absolute side difference is a projected geometry proxy, not clinical prevalence.')
    add('palpebral_length', ['^eye_fissure_length_(r|l|mean)$'], {'HP:0000637': 1, 'HP:0012745': -1},
        note='Fissure length divided by bizygomatic width; population/denominator differences prevent exact clinical equivalence.')
    add('palpebral_height', ['^eye_fissure_height_(r|l|mean)$'], {'HP:0045025': -1}, ['HP:0000508'],
        'Projected opening may reflect gaze, expression or ptosis; short horizontal fissures are not assigned to vertical height.')
    add('palpebral_slant', ['^eye_fissure_slant_(r|l|mean)$'], {'HP:0000582': 1, 'HP:0000494': -1},
        note='Positive is outer canthus above inner canthus. A relative mean shift is not the same as absolute up/down slant.',
        relation='signed_anatomical_angle_proxy')
    add('palpebral_aspect', ['^eye_fissure_aspect_(r|l|mean)$'], {}, ['HP:0045025', 'HP:0000581'],
        'Height/length ratio combines two dimensions; a clinical fissure descriptor does not identify its numerator or denominator effect.')
    add('inner_canthal_spacing', ['^inter_canthal_distance$'], {'HP:0000506': 1, 'HP:0000316': 1, 'HP:0000601': -1},
        note='Inner-canthal/bizygomatic ratio is not an age/population-calibrated telecanthus or bony interorbital measurement.')
    add('pupil_spacing', ['^inter_pupillary_distance$'], {'HP:0000316': 1, 'HP:0000601': -1},
        note='Interpupillary/bizygomatic ratio is a related spacing proxy, not clinical hypertelorism classification.')
    add('outer_canthal_span', ['^outer_canthal_distance$'], {}, ['HP:0000637', 'HP:0012745', 'HP:0000316', 'HP:0000601'],
        'Total outer-canthal span combines spacing and each fissure length; neither trait supplies an unambiguous expected direction.')
    add('canthal_to_pupil_ratio', ['^canthal_to_pupillary_ratio$'], {}, ['HP:0000506', 'HP:0000316'],
        'Clinical spacing traits do not validate a ratio of two different spacing measurements.')
    add('epicanthal_angle', ['^eye_epicanthus_angle_(r|l|mean)$'], {}, ['HP:0000286'],
        'Inner-canthal spoke angle does not establish the presence or severity of an epicanthal skin fold.')
    add('upper_lid_iris_distance', ['^eye_upper_lid_to_iris_(r|l|mean)$'], {'HP:0000508': -1},
        note='Iris-center/upper-lid separation is normalized by face width and depends on gaze; clinical margin-reflex distance is different.')
    add('lower_lid_iris_distance', ['^eye_lower_lid_to_iris_(r|l|mean)$'], {}, ['HP:0000656', 'HP:0007655'],
        'Central lower-lid/iris distance does not quantify eversion of the lateral lower eyelid, retraction or ectropion interchangeably.')
    add('eye_projected_area', ['^eye_area_(r|l|mean)$', '^eye_fissure_fill_(r|l|mean)$'], {},
        ['HP:0008050'], 'Projected lid-contour area/fill is not eyeball volume or an independent clinical eye-size diagnosis.')
    add('iris_offset_or_gaze', ['^iris_offset_', '^gaze_asym_'], {}, ['HP:0000486'],
        'Iris displacement and inter-eye coordinate differences depend on gaze/pose; strabismus supplies no universal signed expected value.')
    add('eye_asymmetry', ['^eye_.*_asym$'], {}, ['HP:0000324'],
        'General facial asymmetry does not imply an asymmetry in this specific ocular ratio.')
    add('nose_length', ['^nose_length$'], {'HP:0003189': 1, 'HP:0003196': -1},
        note='Nasion/subnasale length divided by face width; not clinical absolute length or protrusion.')
    add('nasal_bridge_length', ['^nose_bridge_length$'], {'HP:0033142': 1, 'HP:0003194': -1},
        note='Nasion-to-tip/midface is not an isolated bony bridge length.')
    add('nasal_bridge_width', ['^nose_bridge_width$'], {'HP:0000431': 1, 'HP:0000446': -1}, ['HP:0005280', 'HP:0000426'],
        'Bridge width does not measure bridge depression or prominence.')
    add('nasal_tip_width', ['^nose_tip_width$'], {'HP:0000455': 1, 'HP:0011832': -1}, ['HP:0000414'],
        'Frontal normalized tip width is not three-dimensional bulbosity.')
    add('nasal_base_width', ['^nose_base_width$'], {'HP:0000445': 1, 'HP:0000460': -1},
        note='Alar span/bizygomatic ratio is related to nose width, not an absolute clinical measurement.')
    add('nasal_ratio', ['^nose_tip_to_bridge_ratio$', '^nose_to_face_area_ratio$'], {}, ['HP:0000366'],
        'Composite width/area ratios have no validated one-trait clinical sign mapping.')
    add('nasal_tip_shape', ['^nose_tip_midline_dent$'], {}, ['HP:0000456', 'HP:0004132'],
        'Vertical tip-to-lateral-point difference is not a validated measure of bifidity or tip dimpling.')
    add('nasal_alar_height', ['^nose_ala_height_'], {}, ['HP:0000430'],
        'Alar crest/base vertical separation is not alar tissue thickness; no clinical direction is assumed.')
    add('nasal_profile_angles', ['^nasolabial_angle$', '^nose_tip_elevation$', '^ala_flare_angle$'], {}, ['HP:0000463'],
        'Frontal 2D signed distances/angles do not establish profile anteversion, nares orientation or columella/upper-lip angle.')
    add('columella_length', ['^columella_length$'], {'HP:0002000': -1},
        note='Frontal tip/subnasale distance divided by midface height is only a columella-related proxy.')
    add('columella_position', ['^columella_hang_below_ala$'], {'HP:0009765': 1}, ['HP:0010763', 'HP:0012807'],
        'Subnasale-to-alar-base signed distance differs from direct clinical columella insertion or exposure.')
    add('philtrum_length', ['^philtrum_length$'], {'HP:0000343': 1, 'HP:0000322': -1}, ['HP:0000319'],
        'Subnasale/upper-lip length is divided by midface height; smoothness is not length.')
    add('philtrum_named_width', ['^philtrum_width$'], {}, ['HP:0000289', 'HP:0011829'],
        'Feature spans upper-lip peaks, not philtral columns; broad/narrow clinical philtrum supplies no validated direction.')
    add('mouth_width', ['^mouth_width$'], {'HP:0000154': 1, 'HP:0000160': -1},
        note='Commissure span normalized by face width; expression and clinical threshold differences remain.')
    add('mouth_opening', ['^mouth_opening$'], {}, ['HP:0000194'],
        'One-photo inner-lip gap cannot establish habitual open-mouth posture.')
    add('upper_vermilion', ['^upper_vermilion_height$'], {'HP:0000215': 1, 'HP:0000219': -1},
        note='Vermilion height/face width; skin-lip thickness, eversion and expression are different constructs.')
    add('lower_vermilion', ['^lower_vermilion_height$'], {'HP:0000179': 1, 'HP:0010282': -1},
        note='Vermilion height/face width; not outward eversion of the lip.')
    add('total_vermilion', ['^vermilion_total$'], {},
        ['HP:0012471', 'HP:0000233', 'HP:0000215', 'HP:0000219', 'HP:0000179', 'HP:0010282'],
        note='Sum of upper/lower projected heights; a thin upper lip alone supplies no direction for the sum.')
    add('cupid_bow', ['^cupid_bow_drop$'], {}, ['HP:0002263'],
        'Midline-to-peak vertical difference has no calibrated mapping to a clinically accentuated Cupid bow.')
    add('mouth_corner', ['^mouth_corner_drop$'], {'HP:0002714': 1, 'HP:0010805': -1},
        note='Corner/upper-lip-center vertical difference is a projected relative position, not a calibrated clinical trait threshold.')
    add('mouth_shape', ['^mouth_triangularity$', '^mouth_tenting$', '^mouth_upper_curvature_c2$'], {},
        ['HP:0000207', 'HP:0010804', 'HP:0200096'], 'Curve coefficients and peak/corner ratios are not validated clinical tenting or triangular-mouth indices.')
    add('lip_eversion_named_geometry', ['^upper_lip_eversion$', '^lower_lip_eversion$'], {},
        ['HP:0010803', 'HP:0000232'], 'The feature is a signed vertical border separation, not a direct measurement of outward vermilion eversion.')
    add('lip_midline_asymmetry', ['^lip_midline_x_std$'], {}, ['HP:0009941'],
        'Midline coordinate SD depends on pose/expression and does not identify a fixed clinical mouth asymmetry.')
    add('chin_height', ['^chin_height$'], {'HP:0400000': 1, 'HP:0000331': -1},
        note='Mentolabial/menton height divided by face width; not mandibular depth or prognathism.')
    add('chin_width', ['^chin_width$'], {'HP:0011822': 1}, ['HP:0000307'],
        'Chin width is not the same as jaw size, pointedness or mandibular projection.')
    add('chin_angle', ['^chin_pointedness_angle$'], {'HP:0000307': -1},
        note='Angle is at menton; a smaller angle is sharper. A higher named pointedness_angle means a blunter measured angle.',
        relation='anatomical_angle_proxy')
    add('jaw_width', ['^jaw_bigonial_width$'], {'HP:0012802': 1, 'HP:0012801': -1},
        note='Bigonial/bizygomatic ratio is not mandibular AP size, retrognathia or prognathism.')
    add('forehead_height', ['^forehead_height$'], {'HP:0000348': 1}, ['HP:0011220', 'HP:0002007', 'HP:0000294'],
        'Model trichion/forehead height is normalized by face width; prominence is depth and hairline localization is uncertain.')
    add('forehead_width', ['^forehead_width_upper$', '^forehead_width_mid$'], {'HP:0000337': 1, 'HP:0000341': -1},
        ['HP:0011220'], 'Forehead width/face width is not frontal bossing or a prominent forehead without a width description.')
    add('forehead_taper', ['^forehead_taper_ratio$'], {}, ['HP:0000290'],
        'Upper/mid-forehead width ratio has no validated clinical taper direction.')
    add('face_aspect', ['^face_aspect_ratio$'], {'HP:0000283': 1, 'HP:0000276': -1, 'HP:0000275': -1, 'HP:0011219': 1},
        note='Width/height increases for relatively broader/shorter geometry; not direct cranial shape or absolute face length.')
    add('face_circularity', ['^face_roundness$'], {'HP:0000311': 1},
        note='Projected outline circularity is related to round appearance but is not clinically calibrated.')
    add('face_width_ratios', ['^face_width_uniformity$', '^jaw_to_forehead_ratio$'], {}, ['HP:0000325', 'HP:0000283'],
        'A composite ratio cannot be assigned a direction from one broad face-shape descriptor.')
    add('face_triangularity', ['^face_triangularity$'], {'HP:0000325': 1},
        note='Greater relative forehead-to-chin width contrast is a triangularity proxy, not proof of the entire clinical gestalt.')
    add('face_asymmetry', ['^face_asymmetry$'], {'HP:0000324': 1},
        note='Projected mirror-landmark discrepancy is related to, but not equivalent to, clinical facial asymmetry.')
    add('midface_width', ['^midface_width$'], {}, ['HP:0000309', 'HP:0011800'],
        'Transverse width does not directly measure midface retrusion or hypoplasia in depth.')
    add('cheek_geometry', ['^malar_bulge_', '^cheek_area_'], {}, ['HP:0000293', 'HP:0000272'],
        'Projected cheek/oval indices are not validated measures of tissue fullness or malar AP projection.')
    return rows
