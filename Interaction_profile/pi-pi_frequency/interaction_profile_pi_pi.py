# Import the libraries

import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import os
import argparse

def pi_pi_interaction(pi_system_1, pi_system_2, distance_cutoff=8, cosine_angle_cutoff=0.8):

    """
    This function calculates the pi-pi interaction between aromatic residues.

    Concept: We calculated the sp2/π interactions based on the definition of a recent
literature with small modification of the algorithm for efficiency. First we filtered all the
pairs of sp2/π groups by using a cutoff of 8 ˚A on the distance between the center of mass of
the two groups. Second we calculated cosine angles between the normal vectors of the two
plains and only kept the groups with absolute values of cosine angles larger than 0.8. Third
both the plains defined by each group were raised by 1.5 ˚A and the distance between the
center of mass of the two new plains were calculated. The pairs with the center of mass
distances less than 4 ˚A were selected as forming the sp2/π interactions.
    """

    # Calculate the center of mass of the aromatic residues
    pi_system_COM_1 = pi_system_1.center_of_mass()
    pi_system_COM_2 = pi_system_2.center_of_mass()

    # Calculate the distance between the cation and the center of mass of the aromatic residues
    distance = np.linalg.norm(pi_system_COM_1 - pi_system_COM_2)

    if distance <= distance_cutoff:

        # Calculate the normal vector to the pi system
        pi_system_atom_1 = pi_system_1.select_atoms("name CG")
        pi_system_atom_2 = pi_system_1.select_atoms("name CD1")
        pi_system_atom_3 = pi_system_1.select_atoms("name CD2")

        # Calculate the cross product of the vectors to get the normal vector
        v1 = pi_system_atom_1.positions.ravel() - pi_system_atom_2.positions.ravel()
        v2 = pi_system_atom_1.positions.ravel() - pi_system_atom_3.positions.ravel()

        normal_vector_1 = np.cross(v1, v2)/np.linalg.norm(np.cross(v1, v2))

        # Calculate the normal vector to the pi system
        pi_system_atom_4 = pi_system_2.select_atoms("name CG")
        pi_system_atom_5 = pi_system_2.select_atoms("name CD1")
        pi_system_atom_6 = pi_system_2.select_atoms("name CD2")

        # Calculate the cross product of the vectors to get the normal vector
        v3 = pi_system_atom_4.positions.ravel() - pi_system_atom_5.positions.ravel()
        v4 = pi_system_atom_4.positions.ravel() - pi_system_atom_6.positions.ravel()

        normal_vector_2 = np.cross(v3, v4)/np.linalg.norm(np.cross(v3, v4))

        # Calculate the cosine angle between the normal vector and the vector between the cation and the center of mass of the aromatic residues
        cosine_angle = np.dot(normal_vector_1, normal_vector_2)/(np.linalg.norm(normal_vector_1)*np.linalg.norm(normal_vector_2))

        if cosine_angle >= cosine_angle_cutoff:

            # Raise the plains by 1.5 ˚A
            pi_system_COM_1 += 1.5*normal_vector_1
            pi_system_COM_2 += 1.5*normal_vector_2

            # Calculate the distance between the center of mass of the two new plains
            distance_new = np.linalg.norm(pi_system_COM_1 - pi_system_COM_2)

            if distance_new <= 4:
                return 1
            else:
                return 0
        
def sp2_pi_interaction(sp2_system, pi_system, distance_cutoff=8, cosine_angle_cutoff=0.8):
    
        """
        This function calculates the sp2-pi interaction between aromatic residues.
        """

        # Calculate the center of mass of the aromatic residues
        sp2_system_COM = sp2_system.center_of_mass()
        pi_system_COM = pi_system.center_of_mass()

        # Calculate the distance between the cation and the center of mass of the aromatic residues
        distance = np.linalg.norm(sp2_system_COM - pi_system_COM)

        if distance <= distance_cutoff:

            # Calculate the normal vector to the sp2 system eg. ARG
            sp2_system_atom_1 = sp2_system.select_atoms("name CZ")
            sp2_system_atom_2 = sp2_system.select_atoms("name NH1")
            sp2_system_atom_3 = sp2_system.select_atoms("name NH2")

            # Calculate the normal vector to the pi system
            pi_system_atom_1 = pi_system.select_atoms("name CG")
            pi_system_atom_2 = pi_system.select_atoms("name CD1")
            pi_system_atom_3 = pi_system.select_atoms("name CD2")

            # Calculate the cross product of the vectors to get the normal vector
            v1 = sp2_system_atom_1.positions.ravel() - sp2_system_atom_2.positions.ravel()
            v2 = sp2_system_atom_1.positions.ravel() - sp2_system_atom_3.positions.ravel()

            normal_vector_1 = np.cross(v1, v2)/np.linalg.norm(np.cross(v1, v2))

            # Calculate the cross product of the vectors to get the normal vector

            v3 = pi_system_atom_1.positions.ravel() - pi_system_atom_2.positions.ravel()
            v4 = pi_system_atom_1.positions.ravel() - pi_system_atom_3.positions.ravel()

            normal_vector_2 = np.cross(v3, v4)/np.linalg.norm(np.cross(v3, v4))

            # Calculate the cosine angle between the normal vector and the vector between the cation and the center of mass of the aromatic residues
            cosine_angle = np.dot(normal_vector_1, normal_vector_2)/(np.linalg.norm(normal_vector_1)*np.linalg.norm(normal_vector_2))

            if cosine_angle >= cosine_angle_cutoff:

                # Raise the plains by 1.5 ˚A    
                sp2_system_COM += 1.5*normal_vector_1
                pi_system_COM += 1.5*normal_vector_2

                # Calculate the distance between the center of mass of the two new plains
                distance_new = np.linalg.norm(sp2_system_COM - pi_system_COM)

                if distance_new <= 4:
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

        # Pi-pi interaction between TYR and TYR

        # Select the aromatic residues
        TYR_aromatic_system_1 = u.select_atoms("resname TYR and name CG CD1 CD2 CE1 CE2 CZ")
        TYR_aromatic_system_2 = u.select_atoms("resname TYR and name CG CD1 CD2 CE1 CE2 CZ")

        # Initialize the number of pi-pi interactions
        TYR_TYR_pi_pi = 0
        total_TYR_TYR_bonds = (len(TYR_aromatic_system_1)/6) * ((len(TYR_aromatic_system_1)/6) -1)
        num_frames = u.trajectory.n_frames

        for ts in u.trajectory:
            for i in range(len(TYR_aromatic_system_1)//6):
                j = i*6
                for k in range(len(TYR_aromatic_system_2)//6):
                    l = k*6
                    if i != k:
                        if pi_pi_interaction(TYR_aromatic_system_1[j:j+6], TYR_aromatic_system_2[l:l+6]):
                            TYR_TYR_pi_pi += 1

        # Calculate the pi-pi interaction profile
        TYR_TYR_pi_pi_profile = TYR_TYR_pi_pi/(total_TYR_TYR_bonds*num_frames)

        # Pi-pi interaction between TYR and PHE

        # Select the aromatic residues
        TYR_aromatic_system = u.select_atoms("resname TYR and name CG CD1 CD2 CE1 CE2 CZ")
        PHE_aromatic_system = u.select_atoms("resname PHE and name CG CD1 CD2 CE1 CE2 CZ")

        # Initialize the number of pi-pi interactions
        TYR_PHE_pi_pi = 0
        total_TYR_PHE_bonds = (len(TYR_aromatic_system)/6) * (len(PHE_aromatic_system)/6)
        num_frames = u.trajectory.n_frames

        for ts in u.trajectory:
            for i in range(len(TYR_aromatic_system)//6):
                j = i*6
                for k in range(len(PHE_aromatic_system)//6):
                    l = k*6
                    if pi_pi_interaction(TYR_aromatic_system[j:j+6], PHE_aromatic_system[l:l+6]):
                        TYR_PHE_pi_pi += 1
        # Calculate the pi-pi interaction profile
        TYR_PHE_pi_pi_profile = TYR_PHE_pi_pi/(total_TYR_PHE_bonds*num_frames)

        #Pi-pi interaction between PHE and PHE

        # Select the aromatic residues
        PHE_aromatic_system_1 = u.select_atoms("resname PHE and name CG CD1 CD2 CE1 CE2 CZ")
        PHE_aromatic_system_2 = u.select_atoms("resname PHE and name CG CD1 CD2 CE1 CE2 CZ")

        # Initialize the number of pi-pi interactions
        PHE_PHE_pi_pi = 0
        total_PHE_PHE_bonds = (len(PHE_aromatic_system_1)/6) * ((len(PHE_aromatic_system_1)/6) -1)
        num_frames = u.trajectory.n_frames

        for ts in u.trajectory:
            for i in range(len(PHE_aromatic_system_1)//6):
                j = i*6
                for k in range(len(PHE_aromatic_system_2)//6):
                    l = k*6
                if PHE_aromatic_system_1[j:j+6] != PHE_aromatic_system_2[l:l+6]:
                    if pi_pi_interaction(PHE_aromatic_system_1[j:j+6], PHE_aromatic_system_2[l:l+6]):
                        PHE_PHE_pi_pi += 1

        # Calculate the pi-pi interaction profile
        PHE_PHE_pi_pi_profile = PHE_PHE_pi_pi/(total_PHE_PHE_bonds*num_frames)

        # sp2-pi interaction between TYR and ARG

        # Select the aromatic residues
        TYR_aromatic_system = u.select_atoms("resname TYR and name CG CD1 CD2 CE1 CE2 CZ")
        ARG_aromatic_system = u.select_atoms("resname ARG and name CZ NH1 NH2")

        # Initialize the number of sp2-pi interactions
        TYR_ARG_sp2_pi = 0
        total_TYR_ARG_bonds = (len(TYR_aromatic_system)/6) * (len(ARG_aromatic_system)/3)
        num_frames = u.trajectory.n_frames

        for ts in u.trajectory:
            for i in range(len(TYR_aromatic_system)//6):
                j = i*6
                for k in range(len(ARG_aromatic_system)//3):
                    l = k*3
                    if sp2_pi_interaction(ARG_aromatic_system[l:l+3], TYR_aromatic_system[j:j+6]):
                        TYR_ARG_sp2_pi += 1
        # Calculate the sp2-pi interaction profile
        TYR_ARG_sp2_pi_profile = TYR_ARG_sp2_pi/(total_TYR_ARG_bonds*num_frames)

        # Put the interaction profiles in a list

        interaction_profile_pi_pi = [TYR_TYR_pi_pi_profile, TYR_PHE_pi_pi_profile, PHE_PHE_pi_pi_profile, TYR_ARG_sp2_pi_profile]

        # Save the interaction profile
        np.savetxt("interaction_profiles_pi_pi.dat", interaction_profile_pi_pi)

    else:
        print("The file does not exist", flush=True)

    return interaction_profile_pi_pi

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

