# Import the libraries

import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt


# Define the function to calculate the H-bond interaction profile

# Interaction between Serine and Serine
Mut16_ser_oxygen = u.select_atoms("index 0:220399 and resname SER and name OG").positions
Mut16_ser_hydrogen = u.select_atoms("index 0:220399 and resname SER and name HG").positions
Mut8_ser_oxygen = u.select_atoms("index 220400:228450 and resname SER and name OG").positions
Mut8_ser_hydrogen = u.select_atoms("index 220400:228450 and resname SER and name HG").positions

def Ser_Ser_Hbond_interaction_profile(u, donor, hydrogen, acceptor, cutoff_distance=3.5, cutoff_angle=120):
    """
    This function calculates the H-bond interaction profile between two Serine residues.
    
    """
    # Make serine hydrogen-oxygen vectors to put the angle condition

    Mu16_ser_HO_vector = Mut16_ser_hydrogen - Mut16_ser_oxygen
    Mut8_ser_HO_vector = Mut8_ser_hydrogen - Mut8_ser_oxygen
    Mut16_Mut8_ser_HO_vector = Mut16_ser_hydrogen - Mut8_ser_oxygen
    Mut8_Mut16_ser_HO_vector = Mut8_ser_hydrogen - Mut16_ser_oxygen
   

    # Satisfy the distince and angle condition
    if np.linalg.norm(Mut16_ser_oxygen - Mut8_ser_oxygen) < 3.5:

        return 1

    
