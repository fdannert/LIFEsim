# ---------------------- 
# IMPORT STATEMENTS 
# ---------------------- 
import numpy as np 
import lifesim 

# ---------------------- 
# SETUP YIELDWRAPPER 
# ---------------------- 
 
ywrap = lifesim.ScienceYield( 
    config_path='$config_path', 
    catalog_path='$catalog_path', 
    output_path='$output_path', 
    n_cpu=50, 			# number of CPUs to be used for S/N calculation
    cat_from_ppop=True 		# true if catalog_path points to a .txt file 
) 

# ---------------------- 
# RUN SNR_1h 
# ---------------------- 

# runs the snr_1h method for a range of mirror temperatures  
ywrap.run_sweep_snr( 
    option_name='$option_name',
    option_values=$option_values, 		# defines the aperture sizes for which to run the analysis
    run_name='$run_name'
) 

# combines the normal and the maxsep catalogs 
ywrap.combine_catalog_maxsep(source_name='$run_name') 