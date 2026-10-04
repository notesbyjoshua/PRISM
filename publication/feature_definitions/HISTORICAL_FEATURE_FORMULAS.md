# Historical FaceKit feature formulas

Source revision: `e474628f3cb17a83edad59ac6e1d9e5438753032`

Source file: `src/facekit/core/geometric/extractor.py`

All 125 rows were extracted from the exact historical source revision. The table records the historical Python
expression, supporting local definitions, landmarks detected by deterministic perturbation of the exact code,
normalization, units, bilateral rules, and source line location. It does not use the current 120-feature checkout.

The canonicalization shared by every feature centers landmarks at the midpoint of inner-canthal landmarks 133 and
362, rotates the coordinates to make that axis horizontal, and then applies the feature-specific denominator shown
in the CSV. A small `1e-12` term in the source protects divisions from a zero denominator and is preserved in the
recorded expression.

## Extraction summary

- Features in corrected schema: 125
- Unique definitions: 125
- Exact from source: 125
- Newer-only features: 0

## Definitions

| Feature | Type | Source | Normalization | Landmarks |
|---|---|---|---|---:|
| `eb_thickness_r` | area ratio/index | `_feat_eyebrow` L460-L460 | bizygomatic width squared | 12 |
| `eb_thickness_l` | area ratio/index | `_feat_eyebrow` L461-L461 | bizygomatic width squared | 12 |
| `eb_thickness_mean` | bilateral mean | `_feat_eyebrow` L462-L462 | inherited from directly referenced component features | 22 |
| `eb_thickness_asym` | asymmetry metric | `_feat_eyebrow` L463-L464 | inherited from directly referenced component features | 20 |
| `eb_inter_distance` | normalized distance | `_feat_eyebrow` L468-L468 | bizygomatic width | 16 |
| `eb_arch_r` | derived scalar | `_feat_eyebrow` L477-L477 | local denominator or no additional normalization; see exact expression | 7 |
| `eb_arch_l` | derived scalar | `_feat_eyebrow` L477-L477 | local denominator or no additional normalization; see exact expression | 7 |
| `eb_arch_mean` | bilateral mean | `_feat_eyebrow` L478-L478 | inherited from directly referenced component features | 12 |
| `eb_arch_asym` | asymmetry metric | `_feat_eyebrow` L479-L480 | inherited from directly referenced component features | 12 |
| `eb_position_r` | derived scalar | `_feat_eyebrow` L486-L486 | midface height | 10 |
| `eb_position_l` | derived scalar | `_feat_eyebrow` L486-L486 | midface height | 10 |
| `eb_position_mean` | bilateral mean | `_feat_eyebrow` L487-L487 | inherited from directly referenced component features | 16 |
| `eb_position_asym` | asymmetry metric | `_feat_eyebrow` L488-L489 | inherited from directly referenced component features | 14 |
| `eb_medial_flare_r` | local or normalized ratio/index | `_feat_eyebrow` L503-L503 | local denominator or no additional normalization; see exact expression | 7 |
| `eb_medial_flare_l` | local or normalized ratio/index | `_feat_eyebrow` L503-L503 | local denominator or no additional normalization; see exact expression | 7 |
| `eb_lateral_thickness_r` | area ratio/index | `_feat_eyebrow` L512-L512 | bizygomatic width squared | 12 |
| `eb_lateral_thickness_l` | area ratio/index | `_feat_eyebrow` L512-L512 | bizygomatic width squared | 12 |
| `eb_medial_flare_mean` | bilateral mean | `_feat_eyebrow` L514-L514 | inherited from directly referenced component features | 12 |
| `eb_lateral_thickness_mean` | bilateral mean | `_feat_eyebrow` L515-L516 | inherited from directly referenced component features | 22 |
| `eb_lateral_thickness_asym` | asymmetry metric | `_feat_eyebrow` L517-L519 | inherited from directly referenced component features | 20 |
| `eye_fissure_length_r` | normalized distance | `_feat_eye` L531-L531 | bizygomatic width | 4 |
| `eye_fissure_height_r` | normalized distance | `_feat_eye` L532-L532 | bizygomatic width | 4 |
| `eye_fissure_length_l` | normalized distance | `_feat_eye` L531-L531 | bizygomatic width | 4 |
| `eye_fissure_height_l` | normalized distance | `_feat_eye` L532-L532 | bizygomatic width | 4 |
| `eye_fissure_length_mean` | bilateral mean | `_feat_eye` L533-L533 | inherited from directly referenced component features | 6 |
| `eye_fissure_height_mean` | bilateral mean | `_feat_eye` L534-L534 | inherited from directly referenced component features | 6 |
| `eye_fissure_length_asym` | asymmetry metric | `_feat_eye` L535-L536 | inherited from directly referenced component features | 4 |
| `eye_fissure_height_asym` | asymmetry metric | `_feat_eye` L537-L538 | inherited from directly referenced component features | 4 |
| `eye_fissure_slant_r` | angle | `_feat_eye` L544-L544 | none; angular result in degrees | 3 |
| `eye_fissure_slant_l` | angle | `_feat_eye` L544-L544 | none; angular result in degrees | 3 |
| `eye_fissure_slant_mean` | bilateral mean | `_feat_eye` L545-L545 | inherited from directly referenced component features | 4 |
| `eye_fissure_slant_asym` | asymmetry metric | `_feat_eye` L546-L546 | inherited from directly referenced component features | 4 |
| `eye_fissure_aspect_r` | local or normalized ratio/index | `_feat_eye` L548-L548 | inherited from directly referenced component features | 4 |
| `eye_fissure_aspect_l` | local or normalized ratio/index | `_feat_eye` L549-L549 | inherited from directly referenced component features | 4 |
| `eye_fissure_aspect_mean` | bilateral mean | `_feat_eye` L550-L550 | inherited from directly referenced component features | 8 |
| `inter_canthal_distance` | normalized distance | `_feat_eye` L552-L552 | bizygomatic width | 4 |
| `inter_pupillary_distance` | normalized distance | `_feat_eye` L553-L553 | bizygomatic width | 4 |
| `outer_canthal_distance` | normalized distance | `_feat_eye` L554-L554 | bizygomatic width | 4 |
| `canthal_to_pupillary_ratio` | local or normalized ratio/index | `_feat_eye` L555-L556 | inherited from directly referenced component features | 4 |
| `eye_epicanthus_angle_r` | angle | `_feat_eye` L561-L561 | none; angular result in degrees | 3 |
| `eye_epicanthus_angle_l` | angle | `_feat_eye` L561-L561 | none; angular result in degrees | 3 |
| `eye_epicanthus_angle_mean` | bilateral mean | `_feat_eye` L562-L563 | inherited from directly referenced component features | 6 |
| `eye_epicanthus_angle_asym` | asymmetry metric | `_feat_eye` L564-L566 | inherited from directly referenced component features | 6 |
| `eye_upper_lid_to_iris_r` | derived scalar | `_feat_eye` L572-L572 | local denominator or no additional normalization; see exact expression | 5 |
| `eye_lower_lid_to_iris_r` | derived scalar | `_feat_eye` L573-L573 | local denominator or no additional normalization; see exact expression | 5 |
| `eye_upper_lid_to_iris_l` | derived scalar | `_feat_eye` L572-L572 | local denominator or no additional normalization; see exact expression | 5 |
| `eye_lower_lid_to_iris_l` | derived scalar | `_feat_eye` L573-L573 | local denominator or no additional normalization; see exact expression | 5 |
| `eye_upper_lid_to_iris_mean` | bilateral mean | `_feat_eye` L575-L576 | inherited from directly referenced component features | 8 |
| `eye_upper_lid_to_iris_asym` | asymmetry metric | `_feat_eye` L577-L579 | inherited from directly referenced component features | 8 |
| `eye_lower_lid_to_iris_mean` | bilateral mean | `_feat_eye` L580-L581 | inherited from directly referenced component features | 8 |
| `eye_lower_lid_to_iris_asym` | asymmetry metric | `_feat_eye` L582-L584 | inherited from directly referenced component features | 8 |
| `eye_area_r` | area ratio/index | `_feat_eye` L586-L586 | bizygomatic width squared | 18 |
| `eye_area_l` | area ratio/index | `_feat_eye` L587-L587 | bizygomatic width squared | 18 |
| `eye_area_mean` | bilateral mean | `_feat_eye` L588-L588 | inherited from directly referenced component features | 34 |
| `eye_area_asym` | asymmetry metric | `_feat_eye` L589-L590 | inherited from directly referenced component features | 32 |
| `eye_fissure_fill_r` | local or normalized ratio/index | `_feat_eye` L594-L595 | bizygomatic width squared | 16 |
| `eye_fissure_fill_l` | local or normalized ratio/index | `_feat_eye` L594-L595 | bizygomatic width squared | 16 |
| `eye_fissure_fill_mean` | bilateral mean | `_feat_eye` L596-L596 | inherited from directly referenced component features | 32 |
| `iris_offset_x_r` | local or normalized ratio/index | `_feat_eye` L603-L603 | local denominator or no additional normalization; see exact expression | 4 |
| `iris_offset_y_r` | local or normalized ratio/index | `_feat_eye` L604-L604 | local denominator or no additional normalization; see exact expression | 4 |
| `iris_offset_x_l` | local or normalized ratio/index | `_feat_eye` L603-L603 | local denominator or no additional normalization; see exact expression | 4 |
| `iris_offset_y_l` | local or normalized ratio/index | `_feat_eye` L604-L604 | local denominator or no additional normalization; see exact expression | 4 |
| `gaze_asym_x` | asymmetry metric | `_feat_eye` L605-L605 | inherited from directly referenced component features | 6 |
| `gaze_asym_y` | asymmetry metric | `_feat_eye` L606-L606 | inherited from directly referenced component features | 6 |
| `gaze_asym_norm` | asymmetry metric | `_feat_eye` L607-L607 | inherited from directly referenced component features | 6 |
| `nose_length` | normalized distance | `_feat_nose` L616-L616 | face height | 6 |
| `nose_bridge_length` | normalized distance | `_feat_nose` L617-L617 | midface height | 5 |
| `nose_bridge_width` | normalized distance | `_feat_nose` L618-L618 | bizygomatic width | 4 |
| `nose_tip_width` | normalized distance | `_feat_nose` L619-L619 | bizygomatic width | 4 |
| `nose_base_width` | normalized distance | `_feat_nose` L620-L620 | bizygomatic width | 4 |
| `nose_tip_to_bridge_ratio` | local or normalized ratio/index | `_feat_nose` L622-L622 | inherited from directly referenced component features | 4 |
| `nose_tip_midline_dent` | derived scalar | `_feat_nose` L625-L625 | bizygomatic width | 7 |
| `nose_ala_height_r` | derived scalar | `_feat_nose` L627-L627 | bizygomatic width | 6 |
| `nose_ala_height_l` | derived scalar | `_feat_nose` L628-L628 | bizygomatic width | 6 |
| `nose_ala_height_mean` | bilateral mean | `_feat_nose` L629-L629 | inherited from directly referenced component features | 8 |
| `nose_ala_height_asym` | asymmetry metric | `_feat_nose` L630-L631 | inherited from directly referenced component features | 6 |
| `nostril_region_area` | area ratio/index | `_feat_nose` L635-L635 | bizygomatic width squared | 8 |
| `nasolabial_angle` | angle | `_feat_nose` L637-L637 | none; angular result in degrees | 3 |
| `nose_tip_elevation` | local or normalized ratio/index | `_feat_nose` L638-L638 | bizygomatic width | 6 |
| `ala_flare_angle` | angle | `_feat_nose` L639-L639 | none; angular result in degrees | 3 |
| `columella_length` | derived scalar | `_feat_nose` L642-L642 | midface height | 5 |
| `columella_hang_below_ala` | derived scalar | `_feat_nose` L643-L643 | midface height | 6 |
| `nose_to_face_area_ratio` | area ratio/index | `_feat_nose` L647-L647 | face-oval polygon area | 40 |
| `philtrum_length` | normalized distance | `_feat_philtrum` L655-L655 | midface height | 5 |
| `philtrum_width` | normalized distance | `_feat_philtrum` L656-L656 | mouth width | 4 |
| `mouth_width` | derived scalar | `_feat_mouth_lip` L665-L665 | bizygomatic width | 4 |
| `mouth_opening` | normalized distance | `_feat_mouth_lip` L666-L666 | mouth width | 4 |
| `upper_vermilion_height` | derived scalar | `_feat_mouth_lip` L670-L670 | mouth width | 6 |
| `lower_vermilion_height` | derived scalar | `_feat_mouth_lip` L671-L671 | mouth width | 6 |
| `vermilion_total` | derived scalar | `_feat_mouth_lip` L672-L672 | mouth width | 8 |
| `cupid_bow_drop` | local or normalized ratio/index | `_feat_mouth_lip` L675-L675 | mouth width | 7 |
| `mouth_corner_drop` | local or normalized ratio/index | `_feat_mouth_lip` L678-L678 | mouth width | 5 |
| `mouth_triangularity` | local or normalized ratio/index | `_feat_mouth_lip` L680-L680 | mouth width | 6 |
| `mouth_tenting` | local or normalized ratio/index | `_feat_mouth_lip` L682-L683 | local denominator or no additional normalization; see exact expression | 6 |
| `mouth_upper_curvature_c2` | local or normalized ratio/index | `_feat_mouth_lip` L689-L689 | bizygomatic width | 13 |
| `mouth_upper_fit_residual` | derived scalar | `_feat_mouth_lip` L690-L690 | mouth width | 13 |
| `upper_lip_eversion` | local or normalized ratio/index | `_feat_mouth_lip` L694-L696 | mouth width | 22 |
| `lower_lip_eversion` | local or normalized ratio/index | `_feat_mouth_lip` L699-L701 | mouth width | 22 |
| `lip_outer_to_inner_area` | area ratio/index | `_feat_mouth_lip` L704-L704 | local denominator or no additional normalization; see exact expression | 40 |
| `lip_midline_x_std` | local or normalized ratio/index | `_feat_mouth_lip` L708-L708 | bizygomatic width | 8 |
| `cupid_bow_peak_asym` | asymmetry metric | `_feat_mouth_lip` L709-L710 | mouth width | 6 |
| `chin_height` | derived scalar | `_feat_chin_jaw` L719-L719 | face height | 5 |
| `chin_width` | normalized distance | `_feat_chin_jaw` L720-L720 | bizygomatic width | 4 |
| `chin_pointedness_angle` | angle | `_feat_chin_jaw` L721-L721 | none; angular result in degrees | 3 |
| `jaw_bigonial_width` | normalized distance | `_feat_chin_jaw` L722-L722 | bizygomatic width | 4 |
| `forehead_height` | derived scalar | `_feat_forehead` L731-L731 | face height | 5 |
| `forehead_width_upper` | normalized distance | `_feat_forehead` L732-L732 | bizygomatic width | 4 |
| `forehead_width_mid` | normalized distance | `_feat_forehead` L733-L733 | bizygomatic width | 4 |
| `hairline_height` | derived scalar | `_feat_forehead` L735-L735 | face height | 5 |
| `forehead_taper_ratio` | local or normalized ratio/index | `_feat_forehead` L737-L738 | inherited from directly referenced component features | 4 |
| `face_aspect_ratio` | local or normalized ratio/index | `_feat_face_global` L748-L748 | bizygomatic width | 6 |
| `face_roundness` | area ratio/index | `_feat_face_global` L753-L753 | local denominator or no additional normalization; see exact expression | 36 |
| `face_width_uniformity` | local or normalized ratio/index | `_feat_face_global` L756-L756 | bizygomatic width | 6 |
| `jaw_to_forehead_ratio` | local or normalized ratio/index | `_feat_face_global` L757-L757 | local denominator or no additional normalization; see exact expression | 4 |
| `face_triangularity` | local or normalized ratio/index | `_feat_face_global` L759-L759 | local denominator or no additional normalization; see exact expression | 4 |
| `face_asymmetry` | asymmetry metric | `_feat_face_global` L776-L776 | bizygomatic width | 26 |
| `midface_width` | normalized distance | `_feat_midface_cheek` L785-L785 | bizygomatic width | 4 |
| `malar_bulge_r` | local or normalized ratio/index | `_feat_midface_cheek` L800-L800 | bizygomatic width | 22 |
| `malar_bulge_l` | local or normalized ratio/index | `_feat_midface_cheek` L800-L800 | bizygomatic width | 25 |
| `malar_bulge_mean` | bilateral mean | `_feat_midface_cheek` L802-L802 | inherited from directly referenced component features | 32 |
| `malar_bulge_asym` | asymmetry metric | `_feat_midface_cheek` L803-L804 | inherited from directly referenced component features | 31 |
| `cheek_area_r` | area ratio/index | `_feat_midface_cheek` L812-L812 | bizygomatic width squared | 30 |
| `cheek_area_l` | area ratio/index | `_feat_midface_cheek` L812-L812 | bizygomatic width squared | 28 |
| `cheek_area_mean` | bilateral mean | `_feat_midface_cheek` L813-L813 | inherited from directly referenced component features | 37 |
| `cheek_area_asym` | asymmetry metric | `_feat_midface_cheek` L814-L815 | inherited from directly referenced component features | 37 |

The complete exact expressions and supporting definitions are in `historical_feature_formulas.csv` and
`Supplementary_Table_S6_exact_feature_definitions.csv`.
