# Import the necessary libraries
import MDAnalysis as mda
import mdtraj as md
from MDAnalysis.analysis.distances import distance_array
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.base import AnalysisBase

# Parameters used in the code
chain1 = 100
chain2 = 10
n_res_1 = 140
n_res_2 = 51

# Function to calculate the inter-chain contact map 
def contact_map(uni, chain1, chain2, n_res_1, n_res_2, cutoff=6.0):
    """
    Calculate the contact map between two chain1 and chain2
    in the universe uni.
    Parameters
    ----------
    uni : MDAnalysis Universe
        The universe containing the trajectory
    chain1 : str
        A selection string for the first set of atoms
    chain2 : str
        A selection string for the second set of atoms
    cutoff : float, optional
        The distance cutoff for considering a contact. The default is 6.0.
    Returns
    -------
    contact_map : np.ndarray
        A binary contact map between chain1 and chain2
    """
    # Initialize the contact map as a 2D numpy array
    contact_map = np.zeros((n_res_1, n_res_2), dtype=bool)

    # starting number for the second chains
    start_second_chain = n_res_1 * chain1
    # Write the first for loop to get the residue numbers of the first chain
    for i in range(chain1):
        sel_chain_1 = uni.residues[i*n_res_1:(i+1)*n_res_1]
        # Write the second for loop to get the residue numbers of the second chain
        for j in range(chain2):
            sel_chain_2 = uni.residues[start_second_chain + j*n_res_2:start_second_chain + (j+1)*n_res_2]
            
            # Calculate the pairwise distances between the atoms

            for res_chain_1 in range(n_res_1):
                for res_chain_2 in range(n_res_2):
                    distances = distance_array(sel_chain_1[res_chain_1].atoms.positions, sel_chain_2[res_chain_2].atoms.positions, box=uni.dimensions)
                    
                    # check if the distance between any of the atoms is less than the cutoff

                    in_contact = np.less(distances, cutoff)
                    if np.any(in_contact):
                        contact_map[res_chain_1, res_chain_2] = True

    return contact_map 

# load the trajectory and topology files in MDAnalysis

u = mda.Universe("dynamics.pdb", "dynamics.xtc")

contact_matrix = np.zeros(n_res_1, n_res_2)
for ts in u.trajectory:
    # Calculate the contact map between chain A and chain B
    contact_map_AB = contact_map(u, chain1, chain2, n_res_1, n_res_2, cutoff=6.0)
    contact_matrix = np.sum([contact_matrix, contact_map_AB], axis=0)
contact_matrix = contact_matrix / u.trajectory.n_frames
