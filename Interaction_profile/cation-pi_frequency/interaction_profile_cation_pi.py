# Import the libraries

import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import os
import argparse

def cation_pi_interaction(cation, pi_system, distance_cutoff=6, cosine_angle_cutoff=0.8):
    """
    This function calculates the cation-pi interaction between aromatic residues and cationic residues.

    Concept: Similar to the sp2/π interactions, Cation-π interactions are also
defined by using both a distance and an angle criterion. The distance between the charged
nitrogen in the cationic side chain and the center of mass of the π group is first subjected
to a cutoff of 6 ˚A. The absolute cosine angles between the normal vector of the π plain and
the vector linking the charged nitrogen and the center of mass of the π group is further
subjected to a cutoff of 0.8. The pairs that satisfy both criteria are considered to form the Cation-π interactions.
    """

    # Calculate the center of mass of the aromatic residues
    pi_system_COM = pi_system.center_of_mass()

    # Calculate the distance between the cation and the center of mass of the aromatic residues
    distance = np.linalg.norm(cation - pi_system_COM)

    if distance <= distance_cutoff:

        # Calculate the normal vector to the pi system
        pi_system_atom_1 = pi_system.select_atoms("name CG")
        pi_system_atom_2 = pi_system.select_atoms("name CD1")
        pi_system_atom_3 = pi_system.select_atoms("name CD2")

        # Calculate the cross product of the vectors to get the normal vector
        v1 = pi_system_atom_1.positions - pi_system_atom_2.positions
        v2 = pi_system_atom_1.positions - pi_system_atom_3.positions

        normal_vector = np.cross(v1, v2)/np.linalg.norm(np.cross(v1, v2))

        # Calculate the vector between the cation and the center of mass of the aromatic residues
        vector = cation - pi_system_COM

        # Calculate the cosine angle between the normal vector and the vector between the cation and the center of mass of the aromatic residues
        cosine_angle = np.dot(normal_vector, vector)/(np.linalg.norm(normal_vector)*np.linalg.norm(vector))

        # Subject the distance and the absolute cosine angle to the cutoff
        if cosine_angle >= cosine_angle_cutoff:
            return 1
        else:
            return 0


def main(filename):
    """
    This function calculates the interaction profile between different amino acids.
    
    """
    #Path to the trajectory and topology files
    topology = "/mnt/home/gkumar/Analysis_new2/frame100.gro"
    frame_path = filename

    # Check if the frame file exists
    if os.path.exists(frame_path) and os.path.isfile(frame_path):
    # Create the Universe with the topology and frame files
        u = mda.Universe(topology, frame_path)
        # Now `u` is the MDAnalysis Universe object you can work with
        print(u, flush=True)

        # Cation-pi interaction between ARG and TYR

        # Select the cationic residues
        ARG_cation_1 = u.select_atoms("resname ARG and name NH1")
        ARG_cation_2 = u.select_atoms("resname ARG and name NH2")
        # Select the aromatic residues
        TYR_residue = u.select_atoms("resname TYR and name CG CD1 CD2 CE1 CE2 CZ")

        # Initialize the number of cation-pi interactions
        ARG_TYR_cation_pi = 0
        total_ARG_TYR_bonds = (len(ARG_cation_1)) * (len(TYR_residue)/6)
        num_frames = u.trajectory.n_frames

        # Loop over the trajectory
        for ts in u.trajectory:
            for cation in ARG_cation_1.positions:
                for i in range(len(TYR_residue)//6):
                    j = i*6
                    if cation_pi_interaction(cation, TYR_residue[j:j+6]):
                        ARG_TYR_cation_pi += 1
            for cation in ARG_cation_2.positions:
                for i in range(len(TYR_residue)//6):
                    j = i*6
                    if cation_pi_interaction(cation, TYR_residue[j:j+6]):
                        ARG_TYR_cation_pi += 1
                        
        # Calculate the cation-pi interactions
        ARG_TYR_cation_pi_profile = ARG_TYR_cation_pi/(total_ARG_TYR_bonds*num_frames)

        # Calculate cation-pi interaction between ARG and PHE

        # Select the aromatic residues
        PHE_residue = u.select_atoms("resname PHE and name CG CD1 CD2 CE1 CE2 CZ")

        # Initialize the number of cation-pi interactions
        ARG_PHE_cation_pi = 0
        total_ARG_PHE_bonds = (len(ARG_cation_1)) * (len(PHE_residue)/6)
        num_frames = u.trajectory.n_frames

        # Loop over the trajectory
        for ts in u.trajectory:
            for cation in ARG_cation_1.positions:
                for i in range(len(PHE_residue)//6):
                    j = i*6
                    if cation_pi_interaction(cation, PHE_residue[j:j+6]):
                        ARG_PHE_cation_pi += 1
            for cation in ARG_cation_2.positions:
                for i in range(len(PHE_residue)//6):
                    j = i*6
                    if cation_pi_interaction(cation, PHE_residue[j:j+6]):
                        ARG_PHE_cation_pi += 1

        # Calculate the cation-pi interactions between ARG and PHE
        ARG_PHE_cation_pi_profile = ARG_PHE_cation_pi/(total_ARG_PHE_bonds*num_frames)

        # Calculate cation-pi interaction between LYS and TYR

        # Select the cationic residues
        LYS_cation = u.select_atoms("resname LYS and name NZ")
        # Select the aromatic residues
        TYR_aromatic_system = u.select_atoms("resname TYR and name CG CD1 CD2 CE1 CE2 CZ")

        # Initialize the number of cation-pi interactions
        LYS_TYR_cation_pi = 0
        total_LYS_TYR_bonds = len(LYS_cation) * (len(TYR_aromatic_system)/6)
        num_frames = u.trajectory.n_frames

        # Loop over the trajectory
        for ts in u.trajectory:
            for cation in LYS_cation.positions:
                for i in range(len(TYR_aromatic_system)//6):
                    j = i*6
                    if cation_pi_interaction(cation, TYR_aromatic_system[j:j+6]):
                        LYS_TYR_cation_pi += 1

        # Calculate the cation-pi interactions between LYS and TYR
        LYS_TYR_cation_pi_profile = LYS_TYR_cation_pi/(total_LYS_TYR_bonds*num_frames)

        # Calculate cation-pi interaction between LYS and PHE

        # Select the aromatic residues
        PHE_aromatic_system = u.select_atoms("resname PHE and name CG CD1 CD2 CE1 CE2 CZ")

        # Initialize the number of cation-pi interactions
        LYS_PHE_cation_pi = 0
        total_LYS_PHE_bonds = len(LYS_cation) * (len(PHE_aromatic_system)/6)
        num_frames = u.trajectory.n_frames

        # Loop over the trajectory
        for ts in u.trajectory:
            for cation in LYS_cation.positions:
                for i in range(len(PHE_aromatic_system)//6):
                    j = i*6
                    if cation_pi_interaction(cation, PHE_aromatic_system[j:j+6]):
                        LYS_PHE_cation_pi += 1

        # Calculate the cation-pi interactions between LYS and PHE
        LYS_PHE_cation_pi_profile = LYS_PHE_cation_pi/(total_LYS_PHE_bonds*num_frames)

        # Put the results in a list
        cation_pi_interaction_list = [ARG_TYR_cation_pi_profile, ARG_PHE_cation_pi_profile, LYS_TYR_cation_pi_profile, LYS_PHE_cation_pi_profile]

        #Save the results to a file
        np.savetxt("cation_pi_interaction.dat", cation_pi_interaction_list)

        return cation_pi_interaction_list
    
    else:
        print("The frame file does not exist or is not a file.", flush=True)

if __name__ == "__main__":
    # Create ArgumentParser object
    parser = argparse.ArgumentParser(description="Process a file with MDAnalysis")

    # Add argument for the filename
    parser.add_argument("filename", type=str, help="Path to the file to process")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the main function with the filename
    #print ("calling the function", flush=True)
    main(args.filename)




