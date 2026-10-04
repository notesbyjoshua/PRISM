# Historical feature dependency report

This descriptive map uses the exact FaceKit source revision `e474628f3cb17a83edad59ac6e1d9e5438753032`. It complements the empirical reference-
cohort correlation analysis; it does not alter the historical 6,125-test Benjamini–Hochberg family.

## Landmark use

Across 125 features, landmark counts range from 3 to 40.
The largest dependency sets are:

- `lip_outer_to_inner_area`: 40 landmarks
- `nose_to_face_area_ratio`: 40 landmarks
- `cheek_area_asym`: 37 landmarks
- `cheek_area_mean`: 37 landmarks
- `face_roundness`: 36 landmarks
- `eye_area_mean`: 34 landmarks
- `eye_area_asym`: 32 landmarks
- `eye_fissure_fill_mean`: 32 landmarks
- `malar_bulge_mean`: 32 landmarks
- `malar_bulge_asym`: 31 landmarks

## Shared normalizers

- inherited from directly referenced component features: 39 features
- bizygomatic width: 28 features
- local denominator or no additional normalization; see exact expression: 17 features
- mouth width: 12 features
- bizygomatic width squared: 11 features
- none; angular result in degrees: 7 features
- midface height: 6 features
- face height: 4 features
- face-oval polygon area: 1 features

## Direct derivations

The source contains 93 pairwise edges in which one feature directly references another.

- `eb_thickness_r` ↔ `eb_thickness_mean`
- `eb_thickness_r` ↔ `eb_thickness_asym`
- `eb_thickness_l` ↔ `eb_thickness_mean`
- `eb_thickness_l` ↔ `eb_thickness_asym`
- `eb_thickness_mean` ↔ `eb_thickness_asym`
- `eb_arch_r` ↔ `eb_arch_mean`
- `eb_arch_r` ↔ `eb_arch_asym`
- `eb_arch_l` ↔ `eb_arch_mean`
- `eb_arch_l` ↔ `eb_arch_asym`
- `eb_arch_mean` ↔ `eb_arch_asym`
- `eb_position_r` ↔ `eb_position_mean`
- `eb_position_r` ↔ `eb_position_asym`
- `eb_position_l` ↔ `eb_position_mean`
- `eb_position_l` ↔ `eb_position_asym`
- `eb_position_mean` ↔ `eb_position_asym`
- `eb_medial_flare_r` ↔ `eb_medial_flare_mean`
- `eb_medial_flare_l` ↔ `eb_medial_flare_mean`
- `eb_lateral_thickness_r` ↔ `eb_lateral_thickness_mean`
- `eb_lateral_thickness_r` ↔ `eb_lateral_thickness_asym`
- `eb_lateral_thickness_l` ↔ `eb_lateral_thickness_mean`
- `eb_lateral_thickness_l` ↔ `eb_lateral_thickness_asym`
- `eb_lateral_thickness_mean` ↔ `eb_lateral_thickness_asym`
- `eye_fissure_length_r` ↔ `eye_fissure_length_mean`
- `eye_fissure_length_r` ↔ `eye_fissure_length_asym`
- `eye_fissure_length_r` ↔ `eye_fissure_aspect_r`
- `eye_fissure_height_r` ↔ `eye_fissure_height_mean`
- `eye_fissure_height_r` ↔ `eye_fissure_height_asym`
- `eye_fissure_height_r` ↔ `eye_fissure_aspect_r`
- `eye_fissure_length_l` ↔ `eye_fissure_length_mean`
- `eye_fissure_length_l` ↔ `eye_fissure_length_asym`

## Highest landmark overlap

| Feature 1 | Feature 2 | Shared landmarks |
|---|---|---:|
| `cheek_area_mean` | `cheek_area_asym` | 37 |
| `nose_to_face_area_ratio` | `face_roundness` | 36 |
| `face_roundness` | `cheek_area_asym` | 35 |
| `face_roundness` | `cheek_area_mean` | 35 |
| `nose_to_face_area_ratio` | `cheek_area_asym` | 35 |
| `nose_to_face_area_ratio` | `cheek_area_mean` | 35 |
| `eye_area_asym` | `eye_fissure_fill_mean` | 32 |
| `eye_area_mean` | `eye_area_asym` | 32 |
| `eye_area_mean` | `eye_fissure_fill_mean` | 32 |
| `malar_bulge_mean` | `malar_bulge_asym` | 31 |
| `cheek_area_r` | `cheek_area_asym` | 30 |
| `cheek_area_r` | `cheek_area_mean` | 30 |
| `malar_bulge_mean` | `cheek_area_asym` | 30 |
| `malar_bulge_mean` | `cheek_area_mean` | 30 |
| `face_roundness` | `cheek_area_r` | 29 |
| `malar_bulge_asym` | `cheek_area_asym` | 29 |
| `malar_bulge_asym` | `cheek_area_mean` | 29 |
| `nose_to_face_area_ratio` | `cheek_area_r` | 29 |
| `cheek_area_l` | `cheek_area_asym` | 28 |
| `cheek_area_l` | `cheek_area_mean` | 28 |
| `face_roundness` | `malar_bulge_mean` | 28 |
| `nose_to_face_area_ratio` | `malar_bulge_mean` | 28 |
| `face_roundness` | `cheek_area_l` | 27 |
| `face_roundness` | `malar_bulge_asym` | 27 |
| `nose_to_face_area_ratio` | `cheek_area_l` | 27 |

Complete pairwise edges are in `feature_dependency_edges.csv`; per-feature landmark sets are in
`feature_landmark_dependencies.csv`.
