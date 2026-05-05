yields_path = "/cluster/project/quanz/YOUR_USERNAME/yields"           # path to the yields folder (with subfolders catalogs and runs)
venv_path = "/cluster/home/YOUR_USERNAME/LIFEsim/venv/bin/activate"   # path to the activate script of your virtual environment

today = "20260429"                 # format: YYYYMMDD
catalog_source_date = "20260421"   # format: YYYYMMDD, date of the imported catalogs

option_name =   "primary_temp"        # name of the option to sweep over
option_values = "np.arange(44,50)"    # values to sweep over, as a python expression that can be evaluated with eval()
run_name =      "primarymirror_sweep" # name of the run, used for naming folders and files