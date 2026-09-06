import lifesim
import sys

# Check if an argument was provided
if len(sys.argv) < 2:
    print("Error: No output path provided.")
    sys.exit(1)

# Read output path from the command line argument
# sys.argv[0] is the script name, sys.argv[1] is the first argument
target_output_path = sys.argv[1]

print(f"Starting run with output path: {target_output_path}")

ywrap = lifesim.ScienceYield(
    config_path='$config_path',
    catalog_path='None',
    output_path=target_output_path,  # <--- Using the variable here
    n_cpu=20,
    cat_from_ppop=False
)

ywrap.mange_optimizations(
    scenario_csv='$scenario_csv',
    source_name='ap_merged'
)

ywrap.sweep_mission_time(source_name=None,
                         all_runs=True)
