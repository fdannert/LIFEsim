from pathlib import Path
import shutil
from lifesim.analysis.yield_wrapper import merge_runs

merge_runs(mapping_csv='$mapping_csv',
           merge_csv='$merge_csv',
           output_path='$output_path',)

merge_folder = Path('$output_path2')

catalog_folders_opt = [
    f for f in merge_folder.iterdir()
    if f.is_dir() and f.name not in {"logs", "config_files"}]
catalogs_opt = [(f.name, f.name.split("_", 1)[-1]) for f in catalog_folders_opt]

for full_name, short_name in catalogs_opt:
    catalog_folder = merge_folder / full_name
    output_folder = catalog_folder / "output"
    output_folder.mkdir(exist_ok=True)
    shutil.move(str(catalog_folder / "ap_merged"), str(output_folder / "ap_merged"))

basepath = merge_folder
endpath = "output/"
with open(merge_folder / "config_files" / "optimizer_jobs.csv", "w") as f:
    for full_name, short_name in catalogs_opt:
        f.write(f"{short_name},{basepath}/{full_name}/{endpath}\n")
