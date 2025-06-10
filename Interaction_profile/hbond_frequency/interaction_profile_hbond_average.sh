#!/bin/bash
#SBATCH --job-name=interaction_profile_hbond_average  # Job name
#SBATCH --constraint=icelake
#SBATCH --partition=ccb
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1          # Number of tasks per node
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --time=24:00:00              # Time limit (HH:MM:SS)
#SBATCH --output=interaction_profile_hbond_average.out           # Standard output log file
#SBATCH --error=interaction_profile_hbond_average.err            # Standard error log file

# Load any necessary modules (optional)
source /mnt/home/gkumar/my_env_py/bin/activate


# Run the Python script in the background
python /mnt/home/gkumar/Mut16_Mut8_atomistic_project_analysis/interaction_profile_hbond_average.py 
