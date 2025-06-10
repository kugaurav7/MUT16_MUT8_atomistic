#!/bin/bash
#SBATCH --job-name=interaction_profile_hbond  # Job name
#SBATCH --constraint=icelake
#SBATCH --partition=ccb
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1          # Number of tasks per node
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --time=24:00:00              # Time limit (HH:MM:SS)
#SBATCH --output=interaction_profile_hbond.out           # Standard output log file
#SBATCH --error=interaction_profile_hbond.err            # Standard error log file

# Load any necessary modules (optional)
source /mnt/home/gkumar/my_env_py/bin/activate

# Navigate to the directory where the job should run
#cd /path/to/job/directory

# Run your application or script

# Assign the input argument to the variable i
i=$1

# Initialize count variable
count=0

# Inner loop for CLONE and results directories
for j in {0..99..10}; do
    for k in {0..99}; do
        dir_path="/mnt/ceph/users/gkumar/data_new/RUN$i/CLONE$j/results$k"
        if [ -d "$dir_path" ]; then
            cd "$dir_path" || continue

            # Run the Python script in the background
            python /mnt/home/gkumar/Mut16_Mut8_atomistic_project_analysis/interaction_profile_hbond.py frame"$k".xtc &

            ((count++))

            # Wait for jobs to finish if the count reaches 64
            if [ $count -eq 64 ]; then
                wait
                count=0
            fi
        fi
    done
done

# Wait for any remaining background jobs to finish
wait
