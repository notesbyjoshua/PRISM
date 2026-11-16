"""Targeted source verification; paraphrases, never evidence of literature-wide absence."""
import pandas as pd

ACCESSED = '2026-11-15'
SOURCES = [
    ('NOONAN_2016', 'ORPHA:648', 'Objective studies of the face of Noonan, cardio-facio-cutaneous, and Costello syndromes', 'https://onlinelibrary.wiley.com/doi/full/10.1002/ajmg.a.37736', '2016', 'Noonan age-specific facial description; indexed full-text passage', 'A broad forehead tapering to a small pointed chin is described with age; these remain normalized geometry proxies.'),
    ('RTS_2025', 'ORPHA:783', 'Rubinstein-Taybi syndrome: comprehensive analysis of a Polish cohort', 'https://www.mdpi.com/2073-4425/16/10/1206', '2025', 'Results, 17-person cohort; indexed full-text passage', 'Thick eyebrows are among accompanying craniofacial findings; a source frequency for that individual trait was not extracted.'),
    ('BPTF_2021', 'OMIM:617755', 'Phenotypic expansion of BPTF-related neurodevelopmental disorder with dysmorphic facies and distal limb anomalies', 'https://onlinelibrary.wiley.com/doi/full/10.1002/ajmg.a.62102', '2021', 'Case table; indexed full-text passages', 'Individual cases include broad forehead, triangular face, long nasal bridge and synophrys. These are case-level examples, not a universal gestalt.'),
    ('SMS_GR', 'OMIM:182290', 'GeneReviews: Smith-Magenis syndrome', 'https://www.ncbi.nlm.nih.gov/books/NBK1310/table/sms.T.clinical_features_of_smithmagenis/?report=objectonly', '2026-04-02', 'Table 2; indexed table text', 'Close-spaced eyes are listed among common craniofacial findings; this differs from hypertelorism in the annotation catalogs.'),
    ('MWS_2023', 'OMIM:235730', 'Mutated transcripts of ZEB2 do not undergo nonsense-mediated decay in Mowat-Wilson syndrome', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC10267494/', '2023', 'Case descriptions; indexed full-text passages', 'Two affected individuals have upslanted fissures. This establishes reported variation, not a typical frequency.'),
    ('WBS_2019', 'OMIM:194050', 'Clinical application of chromosomal microarray analysis for Williams-Beuren syndrome in Chinese Han patients', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC6393686/', '2019', 'Results/Table 4; indexed full-text passage', 'Ocular hypertelorism was reported in 24/29 classical cases, alongside short nose in 25/29. This counters treating occasional hypotelorism as universal.'),
    ('WBS_GR', 'OMIM:194050', 'GeneReviews: Williams syndrome', 'https://www.ncbi.nlm.nih.gov/books/NBK1249/', '2023-04-13', 'Suggestive findings; indexed clinical passages', 'Facial descriptors include short nose, broad tip, broad forehead, long philtrum, thick vermilion and wide mouth. Age changes the gestalt.'),
    ('SOTOS_GR', 'OMIM:117550', 'GeneReviews: Sotos syndrome', 'https://www.ncbi.nlm.nih.gov/books/NBK1479/', '2025-06-05', 'Clinical characteristics; indexed clinical passages', 'Tall chin, long narrow face and downslanted fissures are described; chin appearance changes with age.'),
    ('NCBRS_GR', 'OMIM:601358', 'GeneReviews: Nicolaides-Baraitser syndrome', 'https://www.ncbi.nlm.nih.gov/books/NBK321516/', '2025-03-06', 'Clinical characteristics; indexed clinical passages', 'Descriptions include thick eyebrows, downslanted fissures, ptosis, long philtrum, wide mouth, thin upper and thick lower vermilion.'),
    ('KABUKI_GR', 'ORPHA:2322', 'GeneReviews: Kabuki syndrome', 'https://www.ncbi.nlm.nih.gov/books/NBK62111/', '2025-08-14 indexed PDF', 'Facial diagnostic findings; indexed clinical passages', 'Long fissures with lateral lower-lid eversion and broad arched eyebrows with lateral sparsity are clinical findings; central iris gaps are different measurements.'),
    ('CDLS_GR', 'ORPHA:199', 'GeneReviews: Cornelia de Lange syndrome', 'https://www.ncbi.nlm.nih.gov/books/NBK1104/?report=reader', 'Version date not verified', 'Suggestive findings; indexed clinical passages', 'Synophrys, arched thick eyebrows, short nose, long smooth philtrum, thin upper lip and downturned corners are described.'),
    ('BWCFF_GR', 'ORPHA:2995', 'GeneReviews: Baraitser-Winter cerebrofrontofacial syndrome', 'https://www.ncbi.nlm.nih.gov/books/NBK327153/', '2022-03-24', 'Suggestive findings; indexed clinical passages', 'Widely spaced eyes, broad nasal tip, ptosis and arched eyebrows are described for the broad syndrome.'),
    ('JOUBERT_GR', 'ORPHA:140874', 'GeneReviews: Joubert syndrome', 'https://www.ncbi.nlm.nih.gov/books/NBK1325/', '2026-02-12 chapter; indexed PDF may lag', 'Figure 2 caption; indexed clinical passage', 'A pictured child illustrates broad forehead, arched eyebrows and ptosis. Case illustrations do not imply universal traits; neither significant PRISM metric has a direct annotation at the broad parent ID.'),
    ('ZTTK_GR', 'OMIM:617140', 'GeneReviews: ZTTK syndrome', 'https://www.ncbi.nlm.nih.gov/books/NBK618356/', '2025-09-25', 'Clinical description; indexed clinical passages', 'Forehead prominence is described, which does not establish transverse width. Facial asymmetry and palpebral slant vary.'),
    ('ARTHS_2023', 'OMIM:616268', 'Arboleda-Tham syndrome case report', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9859366/', '2023', 'Clinical description; indexed passage', 'Arched eyebrows are reported; lip appearance may change with age. Sparse brows do not imply low geometric arch.'),
    ('ARTHS_2024', 'OMIM:616268', 'Arboleda-Tham syndrome in an Asian girl', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC11109524/', '2024', 'Clinical description; indexed passage', 'Sparse brows, long fissures and hypertelorism broaden the described phenotype; sparsity does not determine brow arch.'),
    ('FBXO11_2018', 'OMIM:618089', 'De novo variants in FBXO11 in 20 individuals with a variable neurodevelopmental disorder', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC6080769/', '2018', 'Original cohort identity; HPO structured frequencies', 'The reviewed catalog attributes downslanted fissures to 5/20. Upward-slant search hits from unrelated genes in other papers were rejected.'),
    ('BCL11A_MEDGEN', 'OMIM:617101', 'NCBI MedGen: Dias-Logan syndrome', 'https://www.ncbi.nlm.nih.gov/medgen/934800', 'Accessed 2026-11-15', 'GeneReviews excerpt and OMIM clinical description; indexed text', 'Facial appearance is variable; downslanted fissures are reported. The structured catalog frequency is 4/8, not a universal descriptor.'),
    ('CLABARS_MEDGEN', 'OMIM:617752', 'NCBI MedGen: Clark-Baraitser syndrome', 'https://www.ncbi.nlm.nih.gov/medgen/443983', 'Accessed 2026-11-15', 'Clinical synopsis; indexed text', 'Facial features vary. Pointed chin in 2/8 and thin upper lip in 3/16 come from structured source annotations, not a universal total-lip-height expectation.'),
    ('FLHS_MEDGEN', 'OMIM:136140', 'NCBI MedGen: Floating-Harbor syndrome', 'https://www.ncbi.nlm.nih.gov/medgen/152667', 'Accessed 2026-11-15', 'Clinical synopsis; indexed text', 'Downturned mouth corners are described, but the relative projected PRISM difference does not diagnose upturned corners.'),
    ('WDSTS_2026', 'OMIM:605130', 'Wiedemann-Steiner syndrome clinical series', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC12862164/', '2026', 'Clinical description; indexed passage', 'Synophrys is reported in a 16-person series; medial landmark spacing is not hair continuity.'),
    ('RTS_2024', 'ORPHA:783', 'Diagnosis and management in Rubinstein-Taybi syndrome: first international consensus statement', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC11137475/', '2024', 'Diagnostic facial criteria; indexed clinical passages', 'The consensus describes nasal/columella gestalt; it does not establish a direction for the FaceKit vertical alar-height metric.'),
    ('FBXW7_2022', 'OMIM:620012', 'Germline variants in FBXW7 are associated with neurodevelopmental disorder', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9069070/', '2022', 'Cohort description; indexed passage', 'Facial findings vary across 35 individuals; no single fixed gestalt establishes directions for every geometric feature.'),
    ('PORRAS_2021', 'OMIM:194050|ORPHA:648', 'Objective differential diagnosis of Noonan and Williams-Beuren syndromes in diverse populations using quantitative facial phenotyping', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8172204/', '2021', 'Full text accessed; Tables 2-3 and methods', 'This compares two syndromes, not syndrome versus unaffected FairFace controls. Its measurement directions cannot be transferred to PRISM without matching the comparator and feature definition.'),
]

# Only explicitly verified positive traits are added. These are not race-specific.
SUPPLEMENTS = [
    ('NOONAN_2016', 'HP:0000337', 'Age-dependent clinical description; frequency not extracted'),
    ('NOONAN_2016', 'HP:0000307', 'Age-dependent clinical description; frequency not extracted'),
    ('RTS_2025', 'HP:0000574', 'Reported in cohort; trait-specific frequency not extracted'),
    ('BPTF_2021', 'HP:0000337', 'Individual case; population frequency not established'),
    ('BPTF_2021', 'HP:0000325', 'Individual case; population frequency not established'),
    ('BPTF_2021', 'HP:0033142', 'Individual case; population frequency not established'),
    ('BPTF_2021', 'HP:0000664', 'Individual case; population frequency not established'),
    ('SMS_GR', 'HP:0000601', 'Common craniofacial finding; >75% table grouping; not an exact metric prevalence'),
    ('MWS_2023', 'HP:0000582', '2 described cases; not a population prevalence'),
    ('WBS_2019', 'HP:0000316', '24/29 classical cases in Chinese Han cohort'),
    ('WBS_2019', 'HP:0003196', '25/29 classical cases in Chinese Han cohort'),
    ('WBS_GR', 'HP:0003196', 'Characteristic descriptor; frequency not extracted'),
    ('WBS_GR', 'HP:0000455', 'Characteristic descriptor; frequency not extracted'),
    ('WBS_GR', 'HP:0000337', 'Characteristic descriptor; frequency not extracted'),
    ('WBS_GR', 'HP:0000343', 'Characteristic descriptor; frequency not extracted'),
    ('WBS_GR', 'HP:0000154', 'Characteristic descriptor; frequency not extracted'),
    ('SOTOS_GR', 'HP:0400000', 'Characteristic age-dependent descriptor; frequency not extracted'),
    ('SOTOS_GR', 'HP:0000276', 'Characteristic descriptor; frequency not extracted'),
    ('SOTOS_GR', 'HP:0000494', 'Characteristic descriptor; frequency not extracted'),
    ('NCBRS_GR', 'HP:0000574', 'Characteristic descriptor; frequency not extracted'),
    ('NCBRS_GR', 'HP:0000494', 'Characteristic descriptor; frequency not extracted'),
    ('NCBRS_GR', 'HP:0000508', 'Characteristic descriptor; frequency not extracted'),
    ('NCBRS_GR', 'HP:0000343', 'Characteristic descriptor; frequency not extracted'),
    ('NCBRS_GR', 'HP:0000154', 'Characteristic descriptor; frequency not extracted'),
    ('NCBRS_GR', 'HP:0000219', 'Characteristic descriptor; frequency not extracted'),
    ('NCBRS_GR', 'HP:0000179', 'Characteristic descriptor; frequency not extracted'),
    ('KABUKI_GR', 'HP:0000637', 'Diagnostic descriptor; frequency not extracted'),
    ('KABUKI_GR', 'HP:0002553', 'Diagnostic descriptor; frequency not extracted'),
    ('KABUKI_GR', 'HP:0005338', 'Diagnostic descriptor; frequency not extracted'),
    ('CDLS_GR', 'HP:0000664', 'Suggestive finding; frequency not extracted'),
    ('CDLS_GR', 'HP:0002553', 'Suggestive finding; frequency not extracted'),
    ('CDLS_GR', 'HP:0000574', 'Suggestive finding; frequency not extracted'),
    ('CDLS_GR', 'HP:0003196', 'Suggestive finding; frequency not extracted'),
    ('CDLS_GR', 'HP:0000343', 'Suggestive finding; frequency not extracted'),
    ('CDLS_GR', 'HP:0000319', 'Suggestive finding; frequency not extracted'),
    ('CDLS_GR', 'HP:0000219', 'Suggestive finding; frequency not extracted'),
    ('CDLS_GR', 'HP:0002714', 'Suggestive finding; frequency not extracted'),
    ('BWCFF_GR', 'HP:0000316', 'Suggestive finding; frequency not extracted'),
    ('BWCFF_GR', 'HP:0000455', 'Suggestive finding; frequency not extracted'),
    ('BWCFF_GR', 'HP:0000508', 'Suggestive finding; frequency not extracted'),
    ('BWCFF_GR', 'HP:0002553', 'Suggestive finding; frequency not extracted'),
]


def source_registry():
    columns = ['source_id', 'disease_ids', 'title', 'url', 'source_version', 'access_and_section', 'review_note']
    frame = pd.DataFrame(SOURCES, columns=columns)
    frame['accessed_utc_date'] = ACCESSED
    frame['literature_absence_can_be_inferred'] = False
    return frame


def supplementary_annotations():
    registry = source_registry().set_index('source_id')
    rows = []
    for source, term, frequency in SUPPLEMENTS:
        item = registry.loc[source]
        rows.append({'disease_id': item.disease_ids, 'disease_name': item.title, 'term_id': term,
                     'frequency': frequency, 'reference': source, 'onset': '', 'sex': '', 'negated': False,
                     'source': 'Targeted_narrative_' + source, 'source_url': item.url})
    return pd.DataFrame(rows)
