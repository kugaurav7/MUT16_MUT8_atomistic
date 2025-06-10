import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import os
import argparse

def salt_bridge_interaction(residue_1, residue_2, distance_cutoff=6):
    """
    This function calculates the salt bridge interaction between two residues.
    Concept: We used a distance cutoff of 6 ˚A on the smallest distance between all charged
nitrogen and oxygen atoms for every pair of charged amino acids to determine the
formation of the salt bridge in our simulations.
    """
    # Calculate the distance between the charged atoms of the residues
    distance = np.linalg.norm(residue_1 - residue_2)
    
    # Subject the distance to the cutoff
    if distance <= distance_cutoff:
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

        # Salt bridge between ARG and ASP
        # Select the cationic residues
        ARG_cation_1 = u.select_atoms("resname ARG and name NH1")
        ARG_cation_2 = u.select_atoms("resname ARG and name NH2")

        # Select the anionic residues
        ASP_anion_1 = u.select_atoms("resname ASP and name OD1")
        ASP_anion_2 = u.select_atoms("resname ASP and name OD2")

        # Initialize the number of salt bridge interactions
        ARG_ASP_salt_bridge = 0
        ARG_ASP_salt_bridge_calc= 0
        total_ARG_ASP_bonds = len(ARG_cation_1) * len(ASP_anion_1)
        num_frames = u.trajectory.n_frames

        # Loop over the trajectory
        for ts in u.trajectory:
            for cation in ARG_cation_1.positions:
                for anion in ASP_anion_1.positions:
                    if salt_bridge_interaction(cation, anion):
                        ARG_ASP_salt_bridge += 1
            
            for cation in ARG_cation_2.positions:
                for anion in ASP_anion_1.positions:
                    if salt_bridge_interaction(cation, anion):
                        ARG_ASP_salt_bridge += 1

            for cation in ARG_cation_2.positions:
                for anion in ASP_anion_2.positions:
                    if salt_bridge_interaction(cation, anion):
                        ARG_ASP_salt_bridge += 1

            for cation in ARG_cation_1.positions:
                for anion in ASP_anion_2.positions:
                    if salt_bridge_interaction(cation, anion):
                        ARG_ASP_salt_bridge += 1

            if ARG_ASP_salt_bridge >= 1:
                ARG_ASP_salt_bridge_calc += 1

        # Calculate the salt bridge interactions
        ARG_ASP_salt_bridge_profile = ARG_ASP_salt_bridge_calc/(total_ARG_ASP_bonds*num_frames)
        
        # Saltbridge between ARG and GLU

        # Select the anionic residues
        GLU_anion_1 = u.select_atoms("resname GLU and name OE1")
        GLU_anion_2 = u.select_atoms("resname GLU and name OE2")

        # Initialize the number of salt bridge interactions
        ARG_GLU_salt_bridge = 0
        ARG_GLU_salt_bridge_calc = 0
        total_ARG_GLU_bonds = len(ARG_cation_1) * len(GLU_anion_1)
        num_frames = u.trajectory.n_frames

        # Loop over the trajectory
        for ts in u.trajectory:
            for cation in ARG_cation_1.positions:
                for anion in GLU_anion_1.positions:
                    if salt_bridge_interaction(cation, anion):
                        ARG_GLU_salt_bridge += 1
            
            for cation in ARG_cation_2.positions:
                for anion in GLU_anion_1.positions:
                    if salt_bridge_interaction(cation, anion):
                        ARG_GLU_salt_bridge += 1

            for cation in ARG_cation_2.positions:
                for anion in GLU_anion_2.positions:
                    if salt_bridge_interaction(cation, anion):
                        ARG_GLU_salt_bridge += 1

            for cation in ARG_cation_1.positions:
                for anion in GLU_anion_2.positions:
                    if salt_bridge_interaction(cation, anion):
                        ARG_GLU_salt_bridge += 1

            if ARG_GLU_salt_bridge >= 1:
                ARG_GLU_salt_bridge_calc += 1

        # Calculate the salt bridge interactions
        ARG_GLU_salt_bridge_profile = ARG_GLU_salt_bridge_calc/(total_ARG_GLU_bonds*num_frames)

        # Saltbridge between LYS and ASP

        # Select the cationic residues
        LYS_cation = u.select_atoms("resname LYS and name NZ")

        # Initialize the number of salt bridge interactions
        LYS_ASP_salt_bridge = 0
        LYS_ASP_salt_bridge_calc = 0
        total_LYS_ASP_bonds = len(ASP_anion_1) * len(LYS_cation)
        num_frames = u.trajectory.n_frames

        # Loop over the trajectory
        for ts in u.trajectory:
            for cation in LYS_cation.positions:
                for anion in ASP_anion_1.positions:
                    if salt_bridge_interaction(cation, anion):
                        LYS_ASP_salt_bridge += 1

            for cation in LYS_cation.positions:
                for anion in ASP_anion_2.positions:
                    if salt_bridge_interaction(cation, anion):
                        LYS_ASP_salt_bridge += 1

            if LYS_ASP_salt_bridge >= 1:
                LYS_ASP_salt_bridge_calc += 1

        # Calculate the salt bridge interactions
        LYS_ASP_salt_bridge_profile = LYS_ASP_salt_bridge_calc/(total_LYS_ASP_bonds*num_frames)

        # Saltbridge between LYS and GLU

        # Initialize the number of salt bridge interactions
        LYS_GLU_salt_bridge = 0
        LYS_GLU_salt_bridge_calc = 0
        total_LYS_GLU_bonds = len(GLU_anion_1) * len(LYS_cation)
        num_frames = u.trajectory.n_frames

        # Loop over the trajectory
        for ts in u.trajectory:
            for cation in LYS_cation.positions:
                for anion in GLU_anion_1.positions:
                    if salt_bridge_interaction(cation, anion):
                        LYS_GLU_salt_bridge += 1

            for cation in LYS_cation.positions:
                for anion in GLU_anion_2.positions:
                    if salt_bridge_interaction(cation, anion):
                        LYS_GLU_salt_bridge += 1

            if LYS_GLU_salt_bridge >= 1:
                LYS_GLU_salt_bridge_calc += 1

        # Calculate the salt bridge interactions
        LYS_GLU_salt_bridge_profile = LYS_GLU_salt_bridge_calc/(total_LYS_GLU_bonds*num_frames)

        # Put the results in a list
        salt_bridge_interaction_list = [ARG_ASP_salt_bridge_profile, ARG_GLU_salt_bridge_profile, LYS_ASP_salt_bridge_profile, LYS_GLU_salt_bridge_profile]

        # Save the results to a file
        np.savetxt("salt_bridge_interaction.dat", salt_bridge_interaction_list)

        return salt_bridge_interaction_list
    
    else:
        print(f"File {frame_path} does not exist", flush=True)
        return None
    
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

