from pathlib import Path
import pandas as pd


# ============================================================
# PATHS — CHANGE IF NEEDED
# ============================================================

INPUT_CSV = (
    "/Users/joshua/Documents/PRISM/data/master_dataset_healthy.csv"
)

OUTPUT_CSV = (
    "/Users/joshua/Documents/PRISM/data/phenotypes_all_healthy.csv"
)


# ============================================================
# LOAD DATA
# ============================================================

df = pd.read_csv(INPUT_CSV)

print(f"Loaded {len(df):,} rows")
print(f"Input columns: {len(df.columns)}")


# ============================================================
# FACEKIT FEATURE COLUMNS
# Exact order requested
# ============================================================

feature_columns = [
    "eb_thickness_r",
    "eb_thickness_l",
    "eb_thickness_mean",
    "eb_thickness_asym",
    "eb_inter_distance",
    "eb_arch_r",
    "eb_arch_l",
    "eb_arch_mean",
    "eb_arch_asym",
    "eb_position_r",
    "eb_position_l",
    "eb_position_mean",
    "eb_position_asym",
    "eb_medial_flare_r",
    "eb_medial_flare_l",
    "eb_lateral_thickness_r",
    "eb_lateral_thickness_l",
    "eb_medial_flare_mean",
    "eb_lateral_thickness_mean",
    "eb_lateral_thickness_asym",

    "eye_fissure_length_r",
    "eye_fissure_height_r",
    "eye_fissure_length_l",
    "eye_fissure_height_l",
    "eye_fissure_length_mean",
    "eye_fissure_height_mean",
    "eye_fissure_length_asym",
    "eye_fissure_height_asym",
    "eye_fissure_slant_r",
    "eye_fissure_slant_l",
    "eye_fissure_slant_mean",
    "eye_fissure_slant_asym",
    "eye_fissure_aspect_r",
    "eye_fissure_aspect_l",
    "eye_fissure_aspect_mean",
    "inter_canthal_distance",
    "inter_pupillary_distance",
    "outer_canthal_distance",
    "canthal_to_pupillary_ratio",
    "eye_epicanthus_angle_r",
    "eye_epicanthus_angle_l",
    "eye_epicanthus_angle_mean",
    "eye_epicanthus_angle_asym",
    "eye_upper_lid_to_iris_r",
    "eye_lower_lid_to_iris_r",
    "eye_upper_lid_to_iris_l",
    "eye_lower_lid_to_iris_l",
    "eye_upper_lid_to_iris_mean",
    "eye_upper_lid_to_iris_asym",
    "eye_lower_lid_to_iris_mean",
    "eye_lower_lid_to_iris_asym",
    "eye_area_r",
    "eye_area_l",
    "eye_area_mean",
    "eye_area_asym",
    "eye_fissure_fill_r",
    "eye_fissure_fill_l",
    "eye_fissure_fill_mean",
    "iris_offset_x_r",
    "iris_offset_y_r",
    "iris_offset_x_l",
    "iris_offset_y_l",
    "gaze_asym_x",
    "gaze_asym_y",
    "gaze_asym_norm",

    "nose_length",
    "nose_bridge_length",
    "nose_bridge_width",
    "nose_tip_width",
    "nose_base_width",
    "nose_tip_to_bridge_ratio",
    "nose_tip_midline_dent",
    "nose_ala_height_r",
    "nose_ala_height_l",
    "nose_ala_height_mean",
    "nose_ala_height_asym",
    "nostril_region_area",
    "nasolabial_angle",
    "nose_tip_elevation",
    "ala_flare_angle",
    "columella_length",
    "columella_hang_below_ala",
    "nose_to_face_area_ratio",

    "philtrum_length",
    "philtrum_width",

    "mouth_width",
    "mouth_opening",
    "upper_vermilion_height",
    "lower_vermilion_height",
    "vermilion_total",
    "cupid_bow_drop",
    "mouth_corner_drop",
    "mouth_triangularity",
    "mouth_tenting",
    "mouth_upper_curvature_c2",
    "mouth_upper_fit_residual",
    "upper_lip_eversion",
    "lower_lip_eversion",
    "lip_outer_to_inner_area",
    "lip_midline_x_std",
    "cupid_bow_peak_asym",

    "chin_height",
    "chin_width",
    "chin_pointedness_angle",
    "jaw_bigonial_width",

    "forehead_height",
    "forehead_width_upper",
    "forehead_width_mid",
    "hairline_height",
    "forehead_taper_ratio",

    "face_aspect_ratio",
    "face_roundness",
    "face_width_uniformity",
    "jaw_to_forehead_ratio",
    "face_triangularity",
    "face_asymmetry",

    "midface_width",
    "malar_bulge_r",
    "malar_bulge_l",
    "malar_bulge_mean",
    "malar_bulge_asym",
    "cheek_area_r",
    "cheek_area_l",
    "cheek_area_mean",
    "cheek_area_asym",
]


# ============================================================
# VERIFY REQUIRED INPUT COLUMNS
# ============================================================

required_input_columns = [
    "file",
    "pose_yaw",
    "pose_pitch",
    "pose_roll",
] + feature_columns

missing = [
    col for col in required_input_columns
    if col not in df.columns
]

if missing:
    print("\nERROR: The following columns are missing:")
    for col in missing:
        print(f"  - {col}")

    raise ValueError(
        f"{len(missing)} required columns are missing from the input CSV."
    )


# ============================================================
# OUTPUT DATASET CREATION
# ============================================================

output = pd.DataFrame(index=df.index)

output["disease"] = "healthy"
output["image_id"] = df["file"].astype(str)
output["frontal_ok"] = True


# ------------------------------------------------------------
# Pose columns
# ------------------------------------------------------------

output["pose_yaw"] = df["pose_yaw"]
output["pose_pitch"] = df["pose_pitch"]
output["pose_roll"] = df["pose_roll"]


# ------------------------------------------------------------
# FaceKit geometric features
# ------------------------------------------------------------

for col in feature_columns:
    output[col] = df[col]


# ============================================================
# SAVE
# ============================================================

output.to_csv(
    OUTPUT_CSV,
    index=False
)


# ============================================================
# CHECK RESULTS
# ============================================================

print()
print("=" * 60)
print("DONE")
print("=" * 60)

print(f"Rows: {len(output):,}")
print(f"Columns: {len(output.columns)}")
print(f"Feature columns: {len(feature_columns)}")

print()
print(f"Saved to:")
print(OUTPUT_CSV)

print()
print("Disease distribution:")
print(output["disease"].value_counts())

print()
print("Frontal distribution:")
print(output["frontal_ok"].value_counts())

print()
print("First 10 columns:")
print(output.columns[:10].tolist())

print()
print("First 5 rows:")
print(
    output[
        [
            "disease",
            "image_id",
            "frontal_ok",
            "pose_yaw",
            "pose_pitch",
            "pose_roll",
        ]
    ].head()
)