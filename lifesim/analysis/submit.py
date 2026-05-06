# imports
import os
import sys
import shutil
import subprocess
from pathlib import Path
from importlib.resources import files
import importlib.util
from string import Template
class BashTemplate(Template):
    delimiter = "@"

def read_template(name: str) -> str:
    """Read a bundled template file from lifesim/analysis/templates/."""
    return (files("lifesim.analysis") / "templates" / name).read_text()

def load_config(config_path: str):
    spec = importlib.util.spec_from_file_location("cfg", config_path)
    mod  = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

# slurm dependencies
def submit(slurm_file, dependency_ids=None):
    """Submit a slurm job and return its job ID. Optionally wait for dependency_ids."""
    cmd = ["sbatch", "--parsable"]
    if dependency_ids:
        dep_str = ":".join(dependency_ids)
        cmd += [f"--dependency=afterok:{dep_str}"]
    cmd.append(str(slurm_file))
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    job_id = result.stdout.strip()
    print(f"  Submitted {slurm_file.name} -> job {job_id}"
          + (f"  (depends on {dependency_ids})" if dependency_ids else ""))
    return job_id

def run(config_path: str):
    """Run the whole workflow with the given configuration file."""
    ###################################
    # general configuration
    ###################################
    cfg = load_config(config_path)

    # personal
    yields              = Path(cfg.yields_path)
    venv_path           = cfg.venv_path

    # run dates
    today               = cfg.today
    catalog_source_date = cfg.catalog_source_date
    queue = cfg.queue

    # sweep options
    option_name         = cfg.option_name
    option_values       = cfg.option_values
    run_name            = cfg.run_name

    # optional configs
    lifesim_config_path   = getattr(cfg, "lifesim_config_path", None)
    optimizer_scenarios_path   = getattr(cfg, "optimizer_scenarios_path", None)
    catalog_merge_path   = getattr(cfg, "catalog_merge_path", None)

    # catalogs
    catalog_folders = [f for f in (yields / "catalogs" / catalog_source_date).iterdir() if f.is_dir()]
    catalogs        = [(f.name, f.name.split("_", 3)[-1]) for f in catalog_folders]

    ###################################
    # snr calculation
    ###################################
    print("\n=== step 1: snr calculation ===")

    template_runyields    = Template(read_template("run_yields_template.py"))
    template_launchscript = Template(read_template("launch_script_template.slurm.sh"))

    for full_name, short_name in catalogs:
        run_folder = yields / "runs" / today / f"{today}_{short_name}"
        run_folder.mkdir(parents=True, exist_ok=True)
        (run_folder / "logs").mkdir(exist_ok=True)
        shutil.rmtree(run_folder / "outputs", ignore_errors=True)
        (run_folder / "output").mkdir(exist_ok=True)

        content = template_runyields.substitute(
            config_path   = yields / "runs" / today / "custom_config.yaml",
            catalog_path  = yields / "catalogs" / catalog_source_date / full_name/ f"{full_name}.txt",
            output_path   = run_folder / "output",
            option_name   = option_name,
            option_values = option_values,
            run_name      = run_name)
        (run_folder / "run_yields.py").write_text(content)

        content = template_launchscript.substitute(
            job_name    = f"{today}_{short_name}",
            output_path = run_folder / "logs" / "python_%j.log",
            python_run  = run_folder / "run_yields.py",
            venv_path   = venv_path,
            queue       = queue)
        (run_folder / "launch_script.slurm.sh").write_text(content)

    content = Path(lifesim_config_path).read_text() if lifesim_config_path else read_template("config_template.yaml")
    (yields / "runs" / today / "custom_config.yaml").write_text(content)

    snr_job_ids = []
    for full_name, short_name in catalogs:
        run_folder = yields / "runs" / today / f"{today}_{short_name}"
        jid = submit(run_folder / "launch_script.slurm.sh")
        snr_job_ids.append(jid)

    print(f"  SNR jobs submitted: {snr_job_ids}")

    ###################################
    # merge 
    ###################################
    print("\n=== step 2: merging ===")

    merge_folder    = yields / "runs" / f"{today}_merge"
    merge_folder.mkdir(parents=True, exist_ok=True)
    (merge_folder / "logs").mkdir(exist_ok=True)
    (merge_folder / "config_files").mkdir(exist_ok=True)

    template_runmerger = Template(read_template("run_merger_template.py"))
    content = template_runmerger.substitute(
            mapping_csv = merge_folder / "config_files" / "catalog_files.csv",
            merge_csv   = merge_folder / "config_files" / "catalog_merge.csv",
            output_path = str(merge_folder) + "/",
            output_path2 = str(merge_folder))
    (merge_folder / "config_files" / "run_merger.py").write_text(content)

    basepath = (yields / "runs" / today)
    endpath = "output/" + str(run_name) + "/"
    with open(merge_folder / "config_files" / "catalog_files.csv", "w") as f:
        f.write("Catalog Name,Catalog Path\n")
        for full_name, short_name in catalogs:
            f.write(f"{short_name},{basepath}/{today}_{short_name}/{endpath}\n")

    if catalog_merge_path is not None:
        raw = Path(catalog_merge_path).read_text()
    else:
        raw = read_template("catalog_merge_template.csv")
    content = Template(raw).safe_substitute(today=today)
    (merge_folder / "config_files" / "catalog_merge.csv").write_text(content, encoding="us-ascii")

    template_launchmerger = Template(read_template("launch_merger_template.slurm.sh"))
    content = template_launchmerger.substitute(
        job_name    = f"{today}_merger",
        output_path = merge_folder / "logs" / "%x_%j.log",
        python_run  = merge_folder / "config_files" / "run_merger.py",
        venv_path   = venv_path,
        queue       = queue)
    (merge_folder / "config_files" / "launch_merger.slurm.sh").write_text(content)

    merge_job_id = submit(
        merge_folder / "config_files" / "launch_merger.slurm.sh",
        dependency_ids=snr_job_ids)

    print(f"  Merge job submitted: {merge_job_id}")

    ###################################
    # optimisation 
    ###################################
    print("\n=== step 3: optimisation ===")

    optimizer_scenarios_destination = merge_folder / "config_files" / "optimizer_scenarios.csv"
    default_optimizer_scenarios_path = files("lifesim.analysis") / "templates" / "optimizer_scenarios_template.csv"
    shutil.copy(optimizer_scenarios_path if optimizer_scenarios_path is not None else default_optimizer_scenarios_path, optimizer_scenarios_destination)

    content = Path(lifesim_config_path).read_text() if lifesim_config_path else read_template("config_template.yaml")
    (merge_folder / "config_files" / "optimizer_config.yaml").write_text(content)

    template_launchopt = Template(read_template("launch_optimizer_template.py"))
    content = template_launchopt.substitute(
            config_path  = merge_folder / "config_files" / "optimizer_config.yaml",
            scenario_csv = merge_folder / "config_files" / "optimizer_scenarios.csv")
    (merge_folder / "config_files" / "launch_optimizer.py").write_text(content)

    template_masterlaunch = BashTemplate(read_template("master_launch_template.slurm.sh"))
    content = template_masterlaunch.substitute(
        job_name     = f"{today}_optimizer",
        output_path  = merge_folder / "logs" / "%x_%A_%a.out",
        optjobs_path = merge_folder / "config_files" / "optimizer_jobs.csv",
        python_run   = merge_folder / "config_files" / "launch_optimizer.py",
        venv_path    = venv_path,
        queue        = queue)
    (merge_folder / "config_files" / "master_launch.slurm.sh").write_text(content)

    opt_job_id = submit(
        merge_folder / "config_files" / "master_launch.slurm.sh",
        dependency_ids=[merge_job_id])
    print(f"  Optimisation job submitted: {opt_job_id}")


command_templates = {
    "init": ("config_example.py", "yield_config.py", "cluster paths and run settings.", None),
    "lifesim_config": ("config_template.yaml", "lifesim_config.yaml", "lifesim and experiment settings. Make sure to keep the file utf-8!", "utf-8"),
    "optimizer_scenarios": ("optimizer_scenarios_template.csv", "optimizer_scenarios.csv", "optimizer scenarios.", None),
    "catalog_merge": ("catalog_merge_template.csv", "catalog_merge.csv", "catalog merge options.", None)
}

def _copy_template(template_name: str, out_name: str, message: str, encoding: str = None):
    out = Path(out_name)
    if out.exists():
        print(f"!!! {out} already exists, not overwriting.")
        return
    content = (files("lifesim.analysis") / "templates" / template_name).read_text(encoding=encoding)
    out.write_text(content, encoding=encoding)
    print(f" Created {out} — edit it with your {message}")


def cli():
    os.chdir(os.environ.get("PWD", os.getcwd()))
    cmd = sys.argv[1] if len(sys.argv) >= 2 else "init"
    if cmd in command_templates:
        _copy_template(*command_templates[cmd])
    else:
        run(cmd)







