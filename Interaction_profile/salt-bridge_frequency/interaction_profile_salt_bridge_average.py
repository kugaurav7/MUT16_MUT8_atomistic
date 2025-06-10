from genericpath import isfile
import numpy as np
import os

# Load the parameters

n_res_1 = 140 # Number of residues in the first protein MUT16
n_res_2 = 51 # Number of residues in the second protein MUT8

# Path to loop over the data files
input_path = "/mnt/ceph/users/gkumar/data_new"
average_interaction_profile_salt_bridge_path = "/mnt/home/gkumar/Mut16_Mut8_atomistic_project_analysis/interaction_profile_average"

# Calculate average of the average and standerd error of mean interaction profile
interaction_profile_salt_bridge_average_average = np.zeros((4), dtype=float)
interaction_profile_salt_bridge_average_sem = np.zeros((4), dtype=float)

# Initialise an empty list to save the interaction profile of each frame
interaction_profile_salt_bridge_list = []

# Loop over the folder to fetch the data files and sun and average the interaction profile

for i in range(0,99):
    counter = 0 # Counter to keep track of the number of frames and calculate the average
    run_path = os.path.join(input_path, f"RUN{i}")
    if not os.path.isdir(run_path):
        continue

    # Initialise the sum and average contact map
    interaction_profile_salt_bridge_sum = np.zeros((4), dtype=float)
    interaction_profile_salt_bridge_average = np.zeros((4), dtype=float)

    for clone in os.listdir(run_path):
        clone_path = os.path.join(run_path, clone)
        if not os.path.isdir(clone_path):
            continue

        for frame in os.listdir(clone_path):
            frame_path = os.path.join(clone_path, frame)
            if not os.path.isdir(frame_path):
                continue

            if os.path.isfile(os.path.join(frame_path, "salt_bridge_interaction.dat")):
                interaction_profile_salt_bridge = np.loadtxt(os.path.join(frame_path, "salt_bridge_interaction.dat"))
                interaction_profile_salt_bridge_sum += interaction_profile_salt_bridge

                counter += 1

    if counter == 0:
        continue

    interaction_profile_salt_bridge_average = interaction_profile_salt_bridge_sum/counter

    # Save the average interaction profile

    np.savetxt(os.path.join(run_path, "interaction_profile_salt_bridge_average.dat"), interaction_profile_salt_bridge_average)

    interaction_profile_salt_bridge_list.append(interaction_profile_salt_bridge_average)

#Convert the list to numpy array
interaction_profile_salt_bridge_list = np.array(interaction_profile_salt_bridge_list)

# Calculate the average and standard error of mean of the interaction profile
interaction_profile_salt_bridge_average_average = np.mean(interaction_profile_salt_bridge_list, axis=0)
interaction_profile_salt_bridge_average_sem = np.std(interaction_profile_salt_bridge_list, axis=0)/np.sqrt(99)

# Save the average interaction profile

np.savetxt(os.path.join(average_interaction_profile_salt_bridge_path, "interaction_profile_salt_bridge_average_average.dat"), interaction_profile_salt_bridge_average_average)
np.savetxt(os.path.join(average_interaction_profile_salt_bridge_path, "interaction_profile_salt_bridge_average_sem.dat"), interaction_profile_salt_bridge_average_sem)


