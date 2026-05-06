#!/bin/bash
#SBATCH --account=$queue
#SBATCH --job-name=$job_name
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=2G
#SBATCH --time=23:30:00
#SBATCH --output=$output_path

module load python  # Or your preferred version
source $venv_path
python $python_run
