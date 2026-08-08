import zipfile
import tempfile
from pathlib import Path, PurePosixPath

import numpy as np
import pandas as pd

from facekit.core.geometric.extractor import GeometricFeatureExtractor


# ============================================================
# CHANGE THESE PATHS
# ============================================================

ZIP_PATH = "/Users/joshua/Documents/PRISM/fairface-img-margin025-trainval.zip"

TRAIN_LABELS = "/Users/joshua/Documents/PRISM/fairface_label_train.csv"
VAL_LABELS = "/Users/joshua/Documents/PRISM/fairface_label_train.csv"

OUTPUT_CSV = (
    "/Users/joshua/Documents/PRISM/data/fairface_1000_age3to9_frontal_balanced.csv"
)


# ============================================================
# SETTINGS
# ============================================================

TOTAL_IMAGES = 1000
TARGET_AGE = "3-9"
RANDOM_SEED = 67

# FaceKit frontal thresholds
extractor = GeometricFeatureExtractor(
    yaw_thresh=15,
    pitch_thresh=15,
    roll_thresh=10,
)


# ============================================================
# LOAD FAIRFACE LABELS
# ============================================================

print("Loading FairFace metadata...")

train = pd.read_csv(TRAIN_LABELS)
val = pd.read_csv(VAL_LABELS)

labels = pd.concat(
    [train, val],
    ignore_index=True
)

print(f"All FairFace labels: {len(labels):,}")


# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

required_columns = {
    "file",
    "race",
    "age",
}

missing = required_columns - set(labels.columns)

if missing:
    raise ValueError(
        f"Missing columns: {missing}\n"
        f"Columns found: {labels.columns.tolist()}"
    )


# ============================================================
# KEEP ONLY CHILDREN AGE 3-9
# ============================================================

labels["age"] = labels["age"].astype(str).str.strip()

labels = labels[
    labels["age"] == TARGET_AGE
].copy()

print()
print(
    f"FairFace images with age {TARGET_AGE}: "
    f"{len(labels):,}"
)

if len(labels) == 0:
    raise ValueError(
        f"No images found with age label {TARGET_AGE!r}."
    )


# ============================================================
# NORMALIZE FILE PATHS
#
# Example:
# train/123.jpg
# val/456.jpg
# ============================================================

def normalize_file_path(x):
    p = PurePosixPath(
        str(x).replace("\\", "/")
    )

    if len(p.parts) >= 2:
        return "/".join(
            p.parts[-2:]
        )

    return str(p)


labels["zip_key"] = (
    labels["file"]
    .apply(normalize_file_path)
)


# ============================================================
# SHOW AVAILABLE RACE COUNTS FOR AGE 3-9
# ============================================================

race_counts = (
    labels["race"]
    .value_counts()
)

races = sorted(
    labels["race"]
    .dropna()
    .unique()
)

print()
print(
    f"Available race counts for age {TARGET_AGE}:"
)

for race in races:
    print(
        f"  {race}: "
        f"{race_counts.get(race, 0):,}"
    )


# ============================================================
# CALCULATE BALANCED QUOTAS
# ============================================================

base = (
    TOTAL_IMAGES
    // len(races)
)

remainder = (
    TOTAL_IMAGES
    % len(races)
)

quotas = {
    race: base
    for race in races
}

# Give the extra +1 slots to groups
# with the most available candidates
extra_races = (
    race_counts
    .sort_values(
        ascending=False
    )
    .index[:remainder]
)

for race in extra_races:
    quotas[race] += 1


print()
print("Target distribution:")

for race in races:
    print(
        f"  {race}: "
        f"{quotas[race]}"
    )

print(
    f"\nTOTAL: "
    f"{sum(quotas.values())}"
)


# ============================================================
# CHECK WHETHER THERE ARE ENOUGH CANDIDATES
# ============================================================

print()
print(
    "Checking whether each race "
    "has enough age 3-9 candidates..."
)

for race in races:

    available = race_counts.get(
        race,
        0
    )

    target = quotas[race]

    if available < target:

        print(
            f"WARNING: {race} has only "
            f"{available} age-{TARGET_AGE} "
            f"images before frontal filtering, "
            f"but target is {target}."
        )


# ============================================================
# OPEN ZIP
# ============================================================

print()
print("Reading ZIP directory...")

with zipfile.ZipFile(
    ZIP_PATH,
    "r"
) as zf:

    # Map:
    #
    # train/1.jpg
    # ->
    # actual file inside ZIP

    zip_lookup = {}

    for member in zf.namelist():

        if member.endswith("/"):
            continue

        p = PurePosixPath(
            member
        )

        if len(p.parts) >= 2:

            key = "/".join(
                p.parts[-2:]
            )

            zip_lookup[key] = member

    print(
        f"Indexed "
        f"{len(zip_lookup):,} "
        f"ZIP files."
    )


    # ========================================================
    # RANDOM NUMBER GENERATOR
    # ========================================================

    rng = np.random.default_rng(
        RANDOM_SEED
    )

    results = []


    # ========================================================
    # PROCESS ONE RACE AT A TIME
    # ========================================================

    for race in races:

        target = quotas[race]

        candidates = labels[
            labels["race"] == race
        ].copy()

        # Randomize candidate order
        order = rng.permutation(
            len(candidates)
        )

        candidates = (
            candidates
            .iloc[order]
        )

        accepted = 0
        attempted = 0
        non_frontal = 0
        no_face = 0
        missing_zip = 0
        errors = 0


        print()
        print("=" * 60)
        print(
            f"Race: {race}"
        )
        print(
            f"Age: {TARGET_AGE}"
        )
        print(
            f"Target: {target}"
        )
        print("=" * 60)


        for _, row in candidates.iterrows():

            if accepted >= target:
                break

            attempted += 1

            zip_key = row["zip_key"]

            member = (
                zip_lookup
                .get(zip_key)
            )

            if member is None:

                missing_zip += 1
                continue


            suffix = (
                Path(zip_key)
                .suffix
            )

            try:

                # ============================================
                # TEMPORARILY EXTRACT ONE IMAGE
                # ============================================

                with tempfile.TemporaryDirectory() as tmpdir:

                    temp_path = (
                        Path(tmpdir)
                        / f"face{suffix}"
                    )

                    with zf.open(
                        member
                    ) as source:

                        with open(
                            temp_path,
                            "wb"
                        ) as destination:

                            destination.write(
                                source.read()
                            )


                    # ========================================
                    # RUN FACEKIT
                    # ========================================

                    feats = extractor.extract(
                        str(temp_path),
                        frontal_check=True
                    )


                # Image has now been deleted automatically


                # ============================================
                # FACE DETECTION FAILED
                # ============================================

                if feats is None:

                    no_face += 1
                    continue


                # ============================================
                # REJECT NON-FRONTAL FACE
                # ============================================

                if (
                    feats.get(
                        "frontal_ok"
                    )
                    is not True
                ):

                    non_frontal += 1
                    continue


                # ============================================
                # ACCEPT
                # ============================================

                record = (
                    row.to_dict()
                )

                record["source"] = (
                    "FairFace"
                )

                record["diagnosis"] = (
                    "Healthy"
                )

                path = PurePosixPath(
                    zip_key
                )

                record["split"] = (
                    path.parent.name
                )

                record["image_id"] = (
                    path.stem
                )

                record.update(
                    feats
                )

                results.append(
                    record
                )

                accepted += 1


                print(
                    f"\r{race}: "
                    f"{accepted}/{target} accepted "
                    f"| {attempted} tested "
                    f"| {non_frontal} non-frontal",
                    end="",
                    flush=True
                )


                # ============================================
                # CHECKPOINT EVERY 25 ACCEPTED FACES
                # ============================================

                if (
                    len(results)
                    % 25 == 0
                ):

                    pd.DataFrame(
                        results
                    ).to_csv(
                        OUTPUT_CSV,
                        index=False
                    )


            except Exception as e:

                errors += 1

                print(
                    f"\nError with "
                    f"{zip_key}: {e}"
                )


        print()
        print()
        print(
            f"{race} finished:"
        )

        print(
            f"  accepted:    {accepted}"
        )

        print(
            f"  attempted:   {attempted}"
        )

        print(
            f"  non-frontal: {non_frontal}"
        )

        print(
            f"  no face:     {no_face}"
        )

        print(
            f"  missing ZIP: {missing_zip}"
        )

        print(
            f"  errors:      {errors}"
        )


        if accepted < target:

            print(
                f"\nWARNING: Only found "
                f"{accepted}/{target} "
                f"valid frontal "
                f"age-{TARGET_AGE} "
                f"faces for {race}."
            )


# ============================================================
# SAVE FINAL DATASET
# ============================================================

df = pd.DataFrame(
    results
)

df.to_csv(
    OUTPUT_CSV,
    index=False
)


# ============================================================
# SUMMARY
# ============================================================

print()
print("=" * 60)
print("DONE")
print("=" * 60)

print(
    f"Final images: "
    f"{len(df):,}"
)

print(
    f"Age: {TARGET_AGE}"
)

print(
    f"Saved to:\n"
    f"{OUTPUT_CSV}"
)


print()
print(
    "Final race distribution:"
)

print(
    df["race"]
    .value_counts()
    .sort_index()
)


print()
print(
    "Age distribution:"
)

print(
    df["age"]
    .value_counts(
        dropna=False
    )
)


print()
print(
    "Frontal status:"
)

print(
    df["frontal_ok"]
    .value_counts(
        dropna=False
    )
)