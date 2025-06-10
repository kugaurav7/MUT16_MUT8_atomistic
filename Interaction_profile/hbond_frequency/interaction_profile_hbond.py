# Import the libraries

import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import os
import argparse

def Hbond_interaction_profile( donor, hydrogen, acceptor, cutoff_distance=3.5, angle_cutoff=120):
    """
    This function calculates the H-bond interaction profile between two Serine residues.
    
    """
    # Satisfy the distince condition

    dist = np.linalg.norm(donor - acceptor)
    cosine = np.dot(hydrogen -donor, hydrogen - acceptor) / (np.linalg.norm(hydrogen - donor) * np.linalg.norm(hydrogen - acceptor))
    cosine_clip = np.clip(cosine, -1, 1)
    angle = np.arccos(cosine_clip) * 180 / np.pi

    if dist < cutoff_distance and angle > angle_cutoff:
        return 1
    
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

        # Interaction between Serine and Serine
        ser_oxygen = u.select_atoms("resname SER and name OG").positions
        ser_hydrogen = u.select_atoms("resname SER and name HG").positions

        ser_ser_hbond = 0
        total_ser_ser_hbond = len(ser_oxygen) * (len(ser_oxygen) - 1) 
        num_frames = len(u.trajectory)

        for ts in u.trajectory:
            for oh1, h1 in zip(ser_oxygen, ser_hydrogen):
                for oh2, h2 in zip(ser_oxygen, ser_hydrogen):
                    if oh1[0] != oh2[0] and oh1[1] != oh2[1] and oh1[2] != oh2[2]:
                        if Hbond_interaction_profile(oh1, h1, oh2) or Hbond_interaction_profile(oh2, h2, oh1):
                            ser_ser_hbond += 1

        ser_ser_hbond_profile = ser_ser_hbond / (total_ser_ser_hbond * num_frames)

        # Interaction between Thr and Thr
        thr_oxygen = u.select_atoms("resname THR and name OG1").positions
        thr_hydrogen = u.select_atoms("resname THR and name HG1").positions

        thr_thr_hbond = 0
        total_thr_thr_hbond = len(thr_oxygen) * (len(thr_oxygen) - 1) 

        for ts in u.trajectory:
            for oh1, h1 in zip(thr_oxygen, thr_hydrogen):
                for oh2, h2 in zip(thr_oxygen, thr_hydrogen):
                    if oh1[0] != oh2[0] and oh1[1] != oh2[1] and oh1[2] != oh2[2]:
                        if Hbond_interaction_profile(oh1, h1, oh2) or Hbond_interaction_profile(oh2, h2, oh1):
                            thr_thr_hbond += 1

        thr_thr_hbond_profile = thr_thr_hbond / (total_thr_thr_hbond * num_frames)

        # Interaction between Asn and Asn
        asn_oxygen = u.select_atoms("resname ASN and name OD1").positions
        asn_hydrogen_1 = u.select_atoms("resname ASN and name HD21").positions
        asn_hydrogen_2 = u.select_atoms("resname ASN and name HD22").positions
        asn_nitrogen = u.select_atoms("resname ASN and name ND2").positions

        asn_asn_hbond = 0
        total_asn_asn_hbond = len(asn_oxygen) * (len(asn_oxygen) - 1) 

        for ts in u.trajectory:
            for oh1, h1, h2, n1 in zip(asn_oxygen, asn_hydrogen_1, asn_hydrogen_2, asn_nitrogen):
                for oh2, h3, h4, n2 in zip(asn_oxygen, asn_hydrogen_1, asn_hydrogen_2, asn_nitrogen):
                    if oh1[0] != oh2[0] and oh1[1] != oh2[1] and oh1[2] != oh2[2]:
                        if Hbond_interaction_profile(oh1, h3, n2) or Hbond_interaction_profile(oh1, h4, n2) or Hbond_interaction_profile(oh2, h1, n1) or Hbond_interaction_profile(oh2, h2, n1):
                            asn_asn_hbond += 1

        asn_asn_hbond_profile = asn_asn_hbond / (total_asn_asn_hbond * num_frames)

        #Interaction between Gln and Gln
        gln_oxygen = u.select_atoms("resname GLN and name OE1").positions
        gln_hydrogen_1 = u.select_atoms("resname GLN and name HE21").positions
        gln_hydrogen_2 = u.select_atoms("resname GLN and name HE22").positions
        gln_nitrogen = u.select_atoms("resname GLN and name NE2").positions

        gln_gln_hbond = 0
        total_gln_gln_hbond = len(gln_oxygen) * (len(gln_oxygen) - 1) 

        for ts in u.trajectory:
            for oh1, h1, h2, n1 in zip(gln_oxygen, gln_hydrogen_1, gln_hydrogen_2, gln_nitrogen):
                for oh2, h3, h4, n2 in zip(gln_oxygen, gln_hydrogen_1, gln_hydrogen_2, gln_nitrogen):
                    if oh1[0] != oh2[0] and oh1[1] != oh2[1] and oh1[2] != oh2[2]:
                        if Hbond_interaction_profile(oh1, h3, n2) or Hbond_interaction_profile(oh1, h4, n2) or Hbond_interaction_profile(oh2, h1, n1) or Hbond_interaction_profile(oh2, h2, n1):
                            gln_gln_hbond += 1

        gln_gln_hbond_profile = gln_gln_hbond / (total_gln_gln_hbond * num_frames)

        #Interaction between Tyr and Tyr
        tyr_oxygen = u.select_atoms("resname TYR and name OH").positions
        tyr_hydrogen = u.select_atoms("resname TYR and name HH").positions

        tyr_tyr_hbond = 0
        total_tyr_tyr_hbond = len(tyr_oxygen) * (len(tyr_oxygen) - 1) 

        for ts in u.trajectory:
            for oh1, h1 in zip(tyr_oxygen, tyr_hydrogen):
                for oh2, h2 in zip(tyr_oxygen, tyr_hydrogen):
                    if oh1[0] != oh2[0] and oh1[1] != oh2[1] and oh1[2] != oh2[2]:
                        if Hbond_interaction_profile(oh1, h1, oh2) or Hbond_interaction_profile(oh2, h2, oh1):
                            tyr_tyr_hbond += 1

        tyr_tyr_hbond_profile = tyr_tyr_hbond / (total_tyr_tyr_hbond * num_frames)

        #Interaction between Ser and Tyr
        ser_tyr_hbond = 0
        total_ser_tyr_hbond = len(ser_oxygen) * len(tyr_oxygen)

        for ts in u.trajectory:
            for oh1, h1 in zip(ser_oxygen, ser_hydrogen):
                for oh2, h2 in zip(tyr_oxygen, tyr_hydrogen):
                    if Hbond_interaction_profile(oh1, h1, oh2) or Hbond_interaction_profile(oh2, h2, oh1):
                        ser_tyr_hbond += 1

        ser_tyr_hbond_profile = ser_tyr_hbond / (total_ser_tyr_hbond * num_frames)

        #Interaction between Thr and Tyr
        thr_tyr_hbond = 0
        total_thr_tyr_hbond = len(thr_oxygen) * len(tyr_oxygen)

        for ts in u.trajectory:
            for oh1, h1 in zip(thr_oxygen, thr_hydrogen):
                for oh2, h2 in zip(tyr_oxygen, tyr_hydrogen):
                    if Hbond_interaction_profile(oh1, h1, oh2) or Hbond_interaction_profile(oh2, h2, oh1):
                        thr_tyr_hbond += 1

        thr_tyr_hbond_profile = thr_tyr_hbond / (total_thr_tyr_hbond * num_frames)

        #Interaction between Asn and Tyr
        asn_tyr_hbond = 0
        total_asn_tyr_hbond = len(asn_oxygen) * len(tyr_oxygen)

        for ts in u.trajectory:
            for oh1 in asn_oxygen:
                for oh2, h1 in zip(tyr_oxygen, tyr_hydrogen):
                    if Hbond_interaction_profile(oh1, h1, oh2):
                        asn_tyr_hbond += 1

        asn_tyr_hbond_profile = asn_tyr_hbond / (total_asn_tyr_hbond * num_frames)

        #Interaction between Gln and Tyr
        gln_tyr_hbond = 0
        total_gln_tyr_hbond = len(gln_oxygen) * len(tyr_oxygen)

        for ts in u.trajectory:
            for oh1 in gln_oxygen:
                for oh2, h1 in zip(tyr_oxygen, tyr_hydrogen):
                    if Hbond_interaction_profile(oh1, h1, oh2):
                        gln_tyr_hbond += 1

        gln_tyr_hbond_profile = gln_tyr_hbond / (total_gln_tyr_hbond * num_frames)

        #Interaction between Ser and Asn
        ser_asn_hbond = 0
        total_ser_asn_hbond = len(ser_oxygen) * len(asn_oxygen)

        for ts in u.trajectory:
            for oh1, h1 in zip(ser_oxygen, ser_hydrogen):
                for oh2 in asn_oxygen:
                    if Hbond_interaction_profile(oh1, h1, oh2):
                        ser_asn_hbond += 1

        ser_asn_hbond_profile = ser_asn_hbond / (total_ser_asn_hbond * num_frames)

        #Interaction between Ser and Gln
        ser_gln_hbond = 0
        total_ser_gln_hbond = len(ser_oxygen) * len(gln_oxygen)

        for ts in u.trajectory:
            for oh1, h1 in zip(ser_oxygen, ser_hydrogen):
                for oh2 in gln_oxygen:
                    if Hbond_interaction_profile(oh1, h1, oh2):
                        ser_gln_hbond += 1

        ser_gln_hbond_profile = ser_gln_hbond / (total_ser_gln_hbond * num_frames)

        #Interaction between Thr and Asn
        thr_asn_hbond = 0
        total_thr_asn_hbond = len(thr_oxygen) * len(asn_oxygen)

        for ts in u.trajectory:
            for oh1, h1 in zip(thr_oxygen, thr_hydrogen):
                for oh2 in asn_oxygen:
                    if Hbond_interaction_profile(oh1, h1, oh2):
                        thr_asn_hbond += 1

        thr_asn_hbond_profile = thr_asn_hbond / total_thr_asn_hbond * num_frames

        #Interaction between Thr and Gln
        thr_gln_hbond = 0
        total_thr_gln_hbond = len(thr_oxygen) * len(gln_oxygen)

        for ts in u.trajectory:
            for oh1, h1 in zip(thr_oxygen, thr_hydrogen):
                for oh2 in gln_oxygen:
                    if Hbond_interaction_profile(oh1, h1, oh2):
                        thr_gln_hbond += 1

        thr_gln_hbond_profile = thr_gln_hbond / total_thr_gln_hbond * num_frames

        # Put all the interaction profiles in a list

        interaction_profiles = [ser_ser_hbond_profile, thr_thr_hbond_profile, asn_asn_hbond_profile, gln_gln_hbond_profile, tyr_tyr_hbond_profile, ser_tyr_hbond_profile, thr_tyr_hbond_profile, asn_tyr_hbond_profile, gln_tyr_hbond_profile, ser_asn_hbond_profile, ser_gln_hbond_profile, thr_asn_hbond_profile, thr_gln_hbond_profile]

        # Save the interaction profiles in a file

        np.savetxt("interaction_profiles.dat", interaction_profiles)

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

