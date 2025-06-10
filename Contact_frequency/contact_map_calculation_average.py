from genericpath import isfile
import numpy as np
import os

# Load the parameters

n_res_1 = 140 # Number of residues in the first protein MUT16
n_res_2 = 51 # Number of residues in the second protein MUT8

# Path to loop over the data files
input_path = "/mnt/ceph/users/gkumar/data_new"
average_contact_maps_path = "/mnt/home/gkumar/Mut16_Mut8_atomistic_project_analysis/average_contact_maps"

# Calculate average of the average contact map and residue-residue contact map
contact_matrix_bb_bb_average_average = np.zeros((n_res_1, n_res_2), dtype=float)
contact_matrix_bb_sc_average_average = np.zeros((n_res_1, n_res_2), dtype=float)
contact_matrix_sc_bb_average_average = np.zeros((n_res_1, n_res_2), dtype=float)
contact_matrix_sc_sc_average_average = np.zeros((n_res_1, n_res_2), dtype=float)

residue_residue_contact_map_bb_bb_normalised_average_average = np.zeros((20, 20), dtype=float)
residue_residue_contact_map_bb_sc_normalised_average_average = np.zeros((20, 20), dtype=float)
residue_residue_contact_map_sc_bb_normalised_average_average = np.zeros((20, 20), dtype=float)
residue_residue_contact_map_sc_sc_normalised_average_average = np.zeros((20, 20), dtype=float)

residue_residue_contact_map_bb_bb_unnormalised_average_average = np.zeros((20, 20), dtype=float)
residue_residue_contact_map_bb_sc_unnormalised_average_average = np.zeros((20, 20), dtype=float)
residue_residue_contact_map_sc_bb_unnormalised_average_average = np.zeros((20, 20), dtype=float)
residue_residue_contact_map_sc_sc_unnormalised_average_average = np.zeros((20, 20), dtype=float)


# Loop over the folder to fetch the data files and sun and average the contact map and residue-residue contact map

for i in range(0,99):
    counter = 0 # Counter to keep track of the number of frames and calculate the average
    run_path = os.path.join(input_path, f"RUN{i}")
    if not os.path.isdir(run_path):
        continue

    # Initialise the sum and average contact map
    contact_matrix_bb_bb_sum = np.zeros((n_res_1, n_res_2), dtype=float)
    contact_matrix_bb_sc_sum = np.zeros((n_res_1, n_res_2), dtype=float)
    contact_matrix_sc_bb_sum = np.zeros((n_res_1, n_res_2), dtype=float)
    contact_matrix_sc_sc_sum = np.zeros((n_res_1, n_res_2), dtype=float)
    contact_matrix_bb_bb_average = np.zeros((n_res_1, n_res_2), dtype=float)
    contact_matrix_bb_sc_average = np.zeros((n_res_1, n_res_2), dtype=float)
    contact_matrix_sc_bb_average = np.zeros((n_res_1, n_res_2), dtype=float)
    contact_matrix_sc_sc_average = np.zeros((n_res_1, n_res_2), dtype=float)

    # Initialise the sum and average residue-residue contact map for normalised
    residue_residue_contact_map_bb_bb_normalized_sum = np.zeros((20, 20), dtype=float)
    residue_residue_contact_map_bb_sc_normalized_sum = np.zeros((20, 20), dtype=float)
    residue_residue_contact_map_sc_bb_normalized_sum = np.zeros((20, 20), dtype=float)
    residue_residue_contact_map_sc_sc_normalized_sum = np.zeros((20, 20), dtype=float)
    residue_residue_contact_map_bb_bb_normalised_average = np.zeros((20, 20), dtype=float)
    residue_residue_contact_map_bb_sc_normalised_average = np.zeros((20, 20), dtype=float)
    residue_residue_contact_map_sc_bb_normalised_average = np.zeros((20, 20), dtype=float)
    residue_residue_contact_map_sc_sc_normalised_average = np.zeros((20, 20), dtype=float)

    # Initialise the sum and average residue-residue contact map for unnormalised
    residue_residue_contact_map_bb_bb_unnormalised_sum = np.zeros((20, 20), dtype=float)
    residue_residue_contact_map_bb_sc_unnormalised_sum = np.zeros((20, 20), dtype=float)
    residue_residue_contact_map_sc_bb_unnormalised_sum = np.zeros((20, 20), dtype=float)
    residue_residue_contact_map_sc_sc_unnormalised_sum = np.zeros((20, 20), dtype=float)
    residue_residue_contact_map_bb_bb_unnormalised_average = np.zeros((20, 20), dtype=float)
    residue_residue_contact_map_bb_sc_unnormalised_average = np.zeros((20, 20), dtype=float)
    residue_residue_contact_map_sc_bb_unnormalised_average = np.zeros((20, 20), dtype=float)
    residue_residue_contact_map_sc_sc_unnormalised_average = np.zeros((20, 20), dtype=float)

    for clone in os.listdir(run_path):
        clone_path = os.path.join(run_path, clone)
        if not os.path.isdir(clone_path):
            continue

        for frame in os.listdir(clone_path):
            frame_path = os.path.join(clone_path, frame)
            if not os.path.isdir(frame_path):
                continue

            if os.path.isfile(os.path.join(frame_path, "contact_matrix_bb_bb.dat")):
                contact_map_bb_bb = np.loadtxt(os.path.join(frame_path, "contact_matrix_bb_bb.dat"))
                contact_map_bb_sc = np.loadtxt(os.path.join(frame_path, "contact_matrix_bb_sc.dat"))
                contact_map_sc_bb = np.loadtxt(os.path.join(frame_path, "contact_matrix_sc_bb.dat"))
                contact_map_sc_sc = np.loadtxt(os.path.join(frame_path, "contact_matrix_sc_sc.dat"))

                residue_residue_contact_map_bb_bb_normalised = np.loadtxt(os.path.join(frame_path, "residue_residue_contact_map_bb_bb_normalised.dat"))
                residue_residue_contact_map_bb_sc_normalised = np.loadtxt(os.path.join(frame_path, "residue_residue_contact_map_bb_sc_normalised.dat"))
                residue_residue_contact_map_sc_bb_normalised = np.loadtxt(os.path.join(frame_path, "residue_residue_contact_map_sc_bb_normalised.dat"))
                residue_residue_contact_map_sc_sc_normalised = np.loadtxt(os.path.join(frame_path, "residue_residue_contact_map_sc_sc_normalised.dat"))
                                                                          
                residue_residue_contact_map_bb_bb_unnormalised = np.loadtxt(os.path.join(frame_path, "residue_residue_contact_map_bb_bb_unnormalised.dat"))
                residue_residue_contact_map_bb_sc_unnormalised = np.loadtxt(os.path.join(frame_path, "residue_residue_contact_map_bb_sc_unnormalised.dat"))
                residue_residue_contact_map_sc_bb_unnormalised = np.loadtxt(os.path.join(frame_path, "residue_residue_contact_map_sc_bb_unnormalised.dat"))
                residue_residue_contact_map_sc_sc_unnormalised = np.loadtxt(os.path.join(frame_path, "residue_residue_contact_map_sc_sc_unnormalised.dat"))

                contact_matrix_bb_bb_sum += contact_map_bb_bb
                contact_matrix_bb_sc_sum += contact_map_bb_sc
                contact_matrix_sc_bb_sum += contact_map_sc_bb
                contact_matrix_sc_sc_sum += contact_map_sc_sc
                
                residue_residue_contact_map_bb_bb_normalized_sum += residue_residue_contact_map_bb_bb_normalised
                residue_residue_contact_map_bb_sc_normalized_sum += residue_residue_contact_map_bb_sc_normalised
                residue_residue_contact_map_sc_bb_normalized_sum += residue_residue_contact_map_sc_bb_normalised
                residue_residue_contact_map_sc_sc_normalized_sum += residue_residue_contact_map_sc_sc_normalised
                
                residue_residue_contact_map_bb_bb_unnormalised_sum += residue_residue_contact_map_bb_bb_unnormalised
                residue_residue_contact_map_bb_sc_unnormalised_sum += residue_residue_contact_map_bb_sc_unnormalised
                residue_residue_contact_map_sc_bb_unnormalised_sum += residue_residue_contact_map_sc_bb_unnormalised
                residue_residue_contact_map_sc_sc_unnormalised_sum += residue_residue_contact_map_sc_sc_unnormalised

                counter += 1
    if counter == 0:
        continue
    contact_matrix_bb_bb_average = contact_matrix_bb_bb_sum / counter
    contact_matrix_bb_sc_average = contact_matrix_bb_sc_sum / counter
    contact_matrix_sc_bb_average = contact_matrix_sc_bb_sum / counter
    contact_matrix_sc_sc_average = contact_matrix_sc_sc_sum / counter

    residue_residue_contact_map_bb_bb_normalised_average = residue_residue_contact_map_bb_bb_normalized_sum / counter
    residue_residue_contact_map_bb_sc_normalised_average = residue_residue_contact_map_bb_sc_normalized_sum / counter
    residue_residue_contact_map_sc_bb_normalised_average = residue_residue_contact_map_sc_bb_normalized_sum / counter
    residue_residue_contact_map_sc_sc_normalised_average = residue_residue_contact_map_sc_sc_normalized_sum / counter

    residue_residue_contact_map_bb_bb_unnormalised_average = residue_residue_contact_map_bb_bb_unnormalised_sum / counter
    residue_residue_contact_map_bb_sc_unnormalised_average = residue_residue_contact_map_bb_sc_unnormalised_sum / counter
    residue_residue_contact_map_sc_bb_unnormalised_average = residue_residue_contact_map_sc_bb_unnormalised_sum / counter
    residue_residue_contact_map_sc_sc_unnormalised_average = residue_residue_contact_map_sc_sc_unnormalised_sum / counter

    # Save the average contact map and residue-residue contact map

    np.savetxt(os.path.join(run_path, "contact_matrix_bb_bb_average.dat") , contact_matrix_bb_bb_average)
    np.savetxt(os.path.join(run_path, "contact_matrix_bb_sc_average.dat") , contact_matrix_bb_sc_average)
    np.savetxt(os.path.join(run_path, "contact_matrix_sc_bb_average.dat") , contact_matrix_sc_bb_average)
    np.savetxt(os.path.join(run_path, "contact_matrix_sc_sc_average.dat") , contact_matrix_sc_sc_average)

    np.savetxt(os.path.join(run_path, "residue_residue_contact_map_bb_bb_normalised_average.dat") , residue_residue_contact_map_bb_bb_normalised_average)
    np.savetxt(os.path.join(run_path, "residue_residue_contact_map_bb_sc_normalised_average.dat") , residue_residue_contact_map_bb_sc_normalised_average)
    np.savetxt(os.path.join(run_path, "residue_residue_contact_map_sc_bb_normalised_average.dat") , residue_residue_contact_map_sc_bb_normalised_average)
    np.savetxt(os.path.join(run_path, "residue_residue_contact_map_sc_sc_normalised_average.dat") , residue_residue_contact_map_sc_sc_normalised_average)

    np.savetxt(os.path.join(run_path, "residue_residue_contact_map_bb_bb_unnormalised_average.dat") , residue_residue_contact_map_bb_bb_unnormalised_average)
    np.savetxt(os.path.join(run_path, "residue_residue_contact_map_bb_sc_unnormalised_average.dat") , residue_residue_contact_map_bb_sc_unnormalised_average)
    np.savetxt(os.path.join(run_path, "residue_residue_contact_map_sc_bb_unnormalised_average.dat") , residue_residue_contact_map_sc_bb_unnormalised_average)
    np.savetxt(os.path.join(run_path, "residue_residue_contact_map_sc_sc_unnormalised_average.dat") , residue_residue_contact_map_sc_sc_unnormalised_average)

    contact_matrix_bb_bb_average_average += contact_matrix_bb_bb_average
    contact_matrix_bb_sc_average_average += contact_matrix_bb_sc_average
    contact_matrix_sc_bb_average_average += contact_matrix_sc_bb_average
    contact_matrix_sc_sc_average_average += contact_matrix_sc_sc_average

    residue_residue_contact_map_bb_bb_normalised_average_average += residue_residue_contact_map_bb_bb_normalised_average
    residue_residue_contact_map_bb_sc_normalised_average_average += residue_residue_contact_map_bb_sc_normalised_average
    residue_residue_contact_map_sc_bb_normalised_average_average += residue_residue_contact_map_sc_bb_normalised_average
    residue_residue_contact_map_sc_sc_normalised_average_average += residue_residue_contact_map_sc_sc_normalised_average

    residue_residue_contact_map_bb_bb_unnormalised_average_average += residue_residue_contact_map_bb_bb_unnormalised_average
    residue_residue_contact_map_bb_sc_unnormalised_average_average += residue_residue_contact_map_bb_sc_unnormalised_average
    residue_residue_contact_map_sc_bb_unnormalised_average_average += residue_residue_contact_map_sc_bb_unnormalised_average
    residue_residue_contact_map_sc_sc_unnormalised_average_average += residue_residue_contact_map_sc_sc_unnormalised_average

contact_matrix_bb_bb_average_average /= 100
contact_matrix_bb_sc_average_average /= 100
contact_matrix_sc_bb_average_average /= 100
contact_matrix_sc_sc_average_average /= 100

residue_residue_contact_map_bb_bb_normalised_average_average /= 100
residue_residue_contact_map_bb_sc_normalised_average_average /= 100
residue_residue_contact_map_sc_bb_normalised_average_average /= 100
residue_residue_contact_map_sc_sc_normalised_average_average /= 100

residue_residue_contact_map_bb_bb_unnormalised_average_average /= 100
residue_residue_contact_map_bb_sc_unnormalised_average_average /= 100
residue_residue_contact_map_sc_bb_unnormalised_average_average /= 100
residue_residue_contact_map_sc_sc_unnormalised_average_average /= 100

np.savetxt(os.path.join(average_contact_maps_path,"contact_matrix_bb_bb_average_average.dat") , contact_matrix_bb_bb_average_average)
np.savetxt(os.path.join(average_contact_maps_path,"contact_matrix_bb_sc_average_average.dat") , contact_matrix_bb_sc_average_average)
np.savetxt(os.path.join(average_contact_maps_path,"contact_matrix_sc_bb_average_average.dat") , contact_matrix_sc_bb_average_average)
np.savetxt(os.path.join(average_contact_maps_path,"contact_matrix_sc_sc_average_average.dat") , contact_matrix_sc_sc_average_average)

np.savetxt(os.path.join(average_contact_maps_path,"residue_residue_contact_map_bb_bb_normalised_average_average.dat") , residue_residue_contact_map_bb_bb_normalised_average_average)
np.savetxt(os.path.join(average_contact_maps_path,"residue_residue_contact_map_bb_sc_normalised_average_average.dat") , residue_residue_contact_map_bb_sc_normalised_average_average)
np.savetxt(os.path.join(average_contact_maps_path,"residue_residue_contact_map_sc_bb_normalised_average_average.dat") , residue_residue_contact_map_sc_bb_normalised_average_average)
np.savetxt(os.path.join(average_contact_maps_path,"residue_residue_contact_map_sc_sc_normalised_average_average.dat") , residue_residue_contact_map_sc_sc_normalised_average_average)

np.savetxt(os.path.join(average_contact_maps_path,"residue_residue_contact_map_bb_bb_unnormalised_average_average.dat") , residue_residue_contact_map_bb_bb_unnormalised_average_average)
np.savetxt(os.path.join(average_contact_maps_path,"residue_residue_contact_map_bb_sc_unnormalised_average_average.dat") , residue_residue_contact_map_bb_sc_unnormalised_average_average)
np.savetxt(os.path.join(average_contact_maps_path,"residue_residue_contact_map_sc_bb_unnormalised_average_average.dat") , residue_residue_contact_map_sc_bb_unnormalised_average_average)
np.savetxt(os.path.join(average_contact_maps_path,"residue_residue_contact_map_sc_sc_unnormalised_average_average.dat") , residue_residue_contact_map_sc_sc_unnormalised_average_average)
