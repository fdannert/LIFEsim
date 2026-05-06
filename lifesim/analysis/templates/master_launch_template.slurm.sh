#!/bin/bash
#SBATCH --account=@queue
#SBATCH --job-name=@job_name
#SBATCH --array=1-7%7        # Run tasks 1 through 7, allow 7 to run at once
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --mem-per-cpu=4G
#SBATCH --time=10:00:00
#SBATCH --output=@output_path # Log files: %A is Job ID, %a is Array Index

# 1. Load Environments
module load python
source @venv_path

# 2. Extract parameters from the CSV file based on the Task ID
# We use 'sed' to grab the specific line number matching the array index
LINE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" @optjobs_path )

# 3. Parse the CSV line
# We set the Internal Field Separator (IFS) to comma to split the line
IFS=',' read -r JOB_LABEL TARGET_PATH <<< "$LINE"

echo "Running Array Task: ${SLURM_ARRAY_TASK_ID}"
echo "Job Label: ${JOB_LABEL}"
echo "Target Path: ${TARGET_PATH}"

# 4. Create the directory if it doesn't exist (optional safety step)
mkdir -p "$TARGET_PATH"

# 5. Run the Python script, passing the path as an argument
python @python_run "$TARGET_PATH"
