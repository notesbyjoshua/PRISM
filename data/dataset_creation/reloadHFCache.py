from datasets import Dataset, concatenate_datasets
from pathlib import Path

cache_dir = Path.home() / (
    ".cache/huggingface/datasets/"
    "GestaltGANwglabshare___gmdb_enhanced_v1.1.0_new/"
    "default/0.0.0/"
    "a94fb856b30eaeeea670b9f06db903afc44e29db"
)

files = sorted(cache_dir.glob("*train-*.arrow"))
print("Arrow files:", len(files))

parts = [Dataset.from_file(str(path)) for path in files]
ds = concatenate_datasets(parts)

print("Number of records:", len(ds))
print("Columns:", ds.column_names)
print("Features:", ds.features)

if "image_restored" in ds.column_names:
    print("Image feature:", ds.features["image_restored"])
    print("First image available:", ds[0]["image_restored"] is not None)