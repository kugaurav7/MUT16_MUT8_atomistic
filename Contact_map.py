# Import the necessary libraries
import MDAnalysis as mda
import mdtraj as md
from MDAnalysis.analysis.distances import distance_array
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.base import AnalysisBase

# Parameters used in the code
chain1 = 100 # Number of chains in the first protein MUT16
chain2 = 10 # Number of chains in the second protein MUT8
n_res_1 = 140 # Number of residues in the first protein MUT16
n_res_2 = 51 # Number of residues in the second protein MUT8

#Path to the trajectory and topology files
topology = "../../../../m2_komet331hpc/kugaurav/Mut16_rde2_interaction/Mut16_condensate/Mut16_Rde2_atomistic_151530/Mut16_Rde2_atomistic_151530/dynamics/dynamics.gro"
trajectory = "../../../../m2_komet331hpc/kugaurav/Mut16_rde2_interaction/Mut16_condensate/Mut16_Rde2_atomistic_151530/Mut16_Rde2_atomistic_151530/dynamics/dynamics.xtc"

# Function to calculate the inter-chain contact map between  MUT16 and  MUT8
def contact_map(uni, chain1, chain2, n_res_1, n_res_2, cutoff=6.0):
    """
    Calculate the contact map between two chain1 and chain2
    in the universe uni.
    Abbreviations:
    -------------
    bb: backbone
    sc: side chain
    ----------

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
    # Initialize the contact map as a 2D numpy array of interaction profile  of MUT16 and MUT8
    contact_map_bb_bb = np.zeros((n_res_1, n_res_2), dtype=bool)
    contact_map_bb_sc = np.zeros((n_res_1, n_res_2), dtype=bool)
    contact_map_sc_bb = np.zeros((n_res_1, n_res_2), dtype=bool)
    contact_map_sc_sc = np.zeros((n_res_1, n_res_2), dtype=bool)

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
                    distances_bb_bb = distance_array(sel_chain_1[res_chain_1].atoms.select_atoms("backbone").positions, sel_chain_2[res_chain_2].atoms.select_atoms("backbone").positions, box=uni.dimensions)
                    distances_bb_sc = distance_array(sel_chain_1[res_chain_1].atoms.select_atoms("backbone").positions, sel_chain_2[res_chain_2].atoms.select_atoms("not backbone").positions, box=uni.dimensions)
                    distances_sc_bb = distance_array(sel_chain_1[res_chain_1].atoms.select_atoms("not backbone").positions, sel_chain_2[res_chain_2].atoms.select_atoms("backbone").positions, box=uni.dimensions)
                    distances_sc_sc = distance_array(sel_chain_1[res_chain_1].atoms.select_atoms("not backbone").positions, sel_chain_2[res_chain_2].atoms.select_atoms("not backbone").positions, box=uni.dimensions)
                    # check if the distance between any of the atoms is less than the cutoff

                    in_contact_bb_bb = np.less(distances_bb_bb, cutoff)
                    if np.any(in_contact_bb_bb):
                        contact_map_bb_bb[res_chain_1, res_chain_2] = True

                    in_contact_bb_sc = np.less(distances_bb_sc, cutoff)
                    if np.any(in_contact_bb_sc):
                        contact_map_bb_sc[res_chain_1, res_chain_2] = True

                    in_contact_sc_bb = np.less(distances_sc_bb, cutoff)
                    if np.any(in_contact_sc_bb):
                        contact_map_sc_bb[res_chain_1, res_chain_2] = True

                    in_contact_sc_sc = np.less(distances_sc_sc, cutoff)
                    if np.any(in_contact_sc_sc):
                        contact_map_sc_sc[res_chain_1, res_chain_2] = True

    return contact_map_bb_bb, contact_map_bb_sc, contact_map_sc_bb, contact_map_sc_sc 


# Write a function to calculate residue-residue contact map for bb:bb, bb:sc, sc:bb, sc:sc
def residue_residue_contact_map(uni, chain1, chain2, n_res_1, n_res_2, contact_matrix_bb_bb, contact_matrix_bb_sc, contact_matrix_sc_bb, contact_matrix_sc_sc):
    """
    Calculate the residue-residue contact map for bb:bb, bb:sc, sc:bb, sc:sc

    Logic:
    ------
    
    
    Returns
    -------
    residue_residue_bb_bb : np.ndarray
        A binary contact map between backbone of chain1 and backbone of chain2
    residue_residue_bb_sc : np.ndarray
        A binary contact map between backbone of chain1 and side chain of chain2
    residue_residue_sc_bb : np.ndarray
        A binary contact map between side chain of chain1 and backbone of chain2
    residue_residue_sc_sc : np.ndarray
        A binary contact map between side chain of chain1 and side chain of chain2
    """
    # Get the sequence of the MUT16 and MUT8 proteins using MDAnalysis
    all_residues = list(uni.residues.resnames)
    MUT16_residues = all_residues[:n_res_1]
    MUT8_residues = all_residues[n_res_1*chain1:n_res_1*chain1+n_res_2]

    residue_residue_matrix_bb_bb = []
    residue_residue_matrix_bb_sc = []   
    residue_residue_matrix_sc_bb = []
    residue_residue_matrix_sc_sc = []

    # Write a for loop to get the residue-residue contact map for bb:bb, bb:sc, sc:bb, sc:sc

    for i, aa_i in enumerate(MUT16_residues):
        # Initialize the row
        row_bb_bb = []
        row_bb_sc = []
        row_sc_bb = []
        row_sc_sc = []

        for j, aa_j in enumerate(MUT8_residues):
            value_bb_bb = contact_matrix_bb_bb[i, j]
            value_bb_sc = contact_matrix_bb_sc[i, j]
            value_sc_bb = contact_matrix_sc_bb[i, j]
            value_sc_sc = contact_matrix_sc_sc[i, j]

            # Resname and values in one list
            element_bb_bb = [aa_i, aa_j, value_bb_bb]
            element_bb_sc = [aa_i, aa_j, value_bb_sc]
            element_sc_bb = [aa_i, aa_j, value_sc_bb]
            element_sc_sc = [aa_i, aa_j, value_sc_sc]

            # Append the element to the row
            row_bb_bb.append(element_bb_bb)
            row_bb_sc.append(element_bb_sc)
            row_sc_bb.append(element_sc_bb)
            row_sc_sc.append(element_sc_sc)

        # Append the row to the matrix
        residue_residue_matrix_bb_bb.append(row_bb_bb)
        residue_residue_matrix_bb_sc.append(row_bb_sc)
        residue_residue_matrix_sc_bb.append(row_sc_bb)
        residue_residue_matrix_sc_sc.append(row_sc_sc)
        
    residue_residue_matrix_bb_bb = np.array(residue_residue_matrix_bb_bb)
    residue_residue_matrix_bb_sc = np.array(residue_residue_matrix_bb_sc)
    residue_residue_matrix_sc_bb = np.array(residue_residue_matrix_sc_bb)
    residue_residue_matrix_sc_sc = np.array(residue_residue_matrix_sc_sc)

    # Create a dictionary to map amino acid names to indices in the matrix
    aa_dict = {'ALA': 0, 'ARG': 1, 'ASN': 2, 'ASP': 3, 'CYS': 4, 'GLU': 5, 'GLN': 6, 'GLY': 7,
           'HIS': 8, 'ILE': 9, 'LEU': 10, 'LYS': 11, 'MET': 12, 'PHE': 13, 'PRO': 14, 'SER': 15,
           'THR': 16, 'TRP': 17, 'TYR': 18, 'VAL': 19}
    
    # Make the residue-residue contact matrix for bb:bb
    aa_matrix_bb_bb = np.zeros((20, 20), dtype=float)
    aa_matrix_count_bb_bb = np.zeros((20, 20), dtype=int)

    for row in residue_residue_matrix_bb_bb:
        for element in row:
            aa_i = element[0]
            aa_j = element[1]
            value = element[2]
            i = aa_dict[aa_i]
            j = aa_dict[aa_j]
            aa_matrix_bb_bb[i, j] += float(value)
            aa_matrix_bb_bb[j, i] += float(value)
            aa_matrix_count_bb_bb[i, j] += 1
            aa_matrix_count_bb_bb[j, i] += 1
                
    aa_matrix_bb_bb_normalised = np.divide(aa_matrix_bb_bb, aa_matrix_count_bb_bb, out=np.zeros_like(aa_matrix_bb_bb), where=aa_matrix_count_bb_bb!=0)
    aa_matrix_bb_bb_unnormalised = aa_matrix_bb_bb

    # Make the residue-residue contact matrix for bb:sc
    aa_matrix_bb_sc = np.zeros((20, 20), dtype=float)
    aa_matrix_count_bb_sc = np.zeros((20, 20), dtype=int)

    for row in residue_residue_matrix_bb_sc:
        for element in row:
            aa_i = element[0]
            aa_j = element[1]
            value = element[2]
            i = aa_dict[aa_i]
            j = aa_dict[aa_j]
            aa_matrix_bb_sc[i, j] += float(value)
            aa_matrix_bb_sc[j, i] += float(value)
            aa_matrix_count_bb_sc[i, j] += 1
            aa_matrix_count_bb_sc[j, i] += 1
            
    aa_matrix_bb_sc_normalised = np.divide(aa_matrix_bb_sc, aa_matrix_count_bb_sc, out=np.zeros_like(aa_matrix_bb_sc), where=aa_matrix_count_bb_sc!=0)
    aa_matrix_bb_sc_unnormalised = aa_matrix_bb_sc

    # Make the residue-residue contact matrix for sc:bb
    aa_matrix_sc_bb = np.zeros((20, 20), dtype=float)
    aa_matrix_count_sc_bb = np.zeros((20, 20), dtype=int)

    for row in residue_residue_matrix_sc_bb:
        for element in row:
            aa_i = element[0]
            aa_j = element[1]
            value = element[2]
            i = aa_dict[aa_i]
            j = aa_dict[aa_j]
            aa_matrix_sc_bb[i, j] += float(value)
            aa_matrix_sc_bb[j, i] += float(value)
            aa_matrix_count_sc_bb[i, j] += 1
            aa_matrix_count_sc_bb[j, i] += 1
            
    aa_matrix_sc_bb_normalised = np.divide(aa_matrix_sc_bb, aa_matrix_count_sc_bb, out=np.zeros_like(aa_matrix_sc_bb), where=aa_matrix_count_sc_bb!=0)
    aa_matrix_sc_bb_unnormalised = aa_matrix_sc_bb

    # Make the residue-residue contact matrix for sc:sc
    aa_matrix_sc_sc = np.zeros((20, 20), dtype=float)
    aa_matrix_count_sc_sc = np.zeros((20, 20), dtype=int)

    for row in residue_residue_matrix_sc_sc:
        for element in row:
            aa_i = element[0]
            aa_j = element[1]
            value = element[2]
            i = aa_dict[aa_i]
            j = aa_dict[aa_j]
            aa_matrix_sc_sc[i, j] += float(value)
            aa_matrix_sc_sc[j, i] += float(value)
            aa_matrix_count_sc_sc[i, j] += 1
            aa_matrix_count_sc_sc[j, i] += 1
            
    aa_matrix_sc_sc_normalised = np.divide(aa_matrix_sc_sc, aa_matrix_count_sc_sc, out=np.zeros_like(aa_matrix_sc_sc), where=aa_matrix_count_sc_sc!=0)
    aa_matrix_sc_sc_unnormalised = aa_matrix_sc_sc

    return aa_matrix_bb_bb_normalised, aa_matrix_bb_bb_unnormalised, aa_matrix_bb_sc_normalised, aa_matrix_bb_sc_unnormalised, aa_matrix_sc_bb_normalised, aa_matrix_sc_bb_unnormalised, aa_matrix_sc_sc_normalised, aa_matrix_sc_sc_unnormalised

# load the trajectory and topology files in MDAnalysis
u = mda.Universe(topology, trajectory)


# Initialize the contact matrix bb:bb, bb:sc, sc:bb, sc:sc
contact_matrix_bb_bb = np.zeros((n_res_1, n_res_2))
contact_matrix_bb_sc = np.zeros((n_res_1, n_res_2))
contact_matrix_sc_bb = np.zeros((n_res_1, n_res_2))
contact_matrix_sc_sc = np.zeros((n_res_1, n_res_2))
contact_matrix_sum_bb_bb = np.zeros((n_res_1, n_res_2))
contact_matrix_sum_bb_sc = np.zeros((n_res_1, n_res_2))
contact_matrix_sum_sc_bb = np.zeros((n_res_1, n_res_2))
contact_matrix_sum_sc_sc = np.zeros((n_res_1, n_res_2))

for ts in u.trajectory[:u.trajectory.n_frames]:
    # Calculate the contact map between bb:bb, bb:sc, sc:bb, sc:sc
    contact_map_function = contact_map(u, chain1, chain2, n_res_1, n_res_2, cutoff=6.0)
    contact_map_bb_bb = contact_map_function[0]
    contact_map_bb_sc = contact_map_function[1]
    contact_map_sc_bb = contact_map_function[2]
    contact_map_sc_sc = contact_map_function[3]

    #Sum the contacts over the trajectory
    contact_matrix_sum_bb_bb = np.sum([contact_matrix_sum_bb_bb, contact_map_bb_bb], axis=0)
    contact_matrix_sum_bb_sc = np.sum([contact_matrix_sum_bb_sc, contact_map_bb_sc], axis=0)
    contact_matrix_sum_sc_bb = np.sum([contact_matrix_sum_sc_bb, contact_map_sc_bb], axis=0)
    contact_matrix_sum_sc_sc = np.sum([contact_matrix_sum_sc_sc, contact_map_sc_sc], axis=0)
    
# Average the contacts over the whole trajectory
contact_matrix_bb_bb = contact_matrix_sum_bb_bb / (u.trajectory.n_frames)
contact_matrix_bb_sc = contact_matrix_sum_bb_sc / (u.trajectory.n_frames)
contact_matrix_sc_bb = contact_matrix_sum_sc_bb / (u.trajectory.n_frames)
contact_matrix_sc_sc = contact_matrix_sum_sc_sc / (u.trajectory.n_frames)

# Save th contact matrix to a .dat file
np.savetxt("contact_matrix_bb_bb.dat", contact_matrix_bb_bb)
np.savetxt("contact_matrix_bb_sc.dat", contact_matrix_bb_sc)
np.savetxt("contact_matrix_sc_bb.dat", contact_matrix_sc_bb)
np.savetxt("contact_matrix_sc_sc.dat", contact_matrix_sc_sc)



# Calculate the residue-residue contact map for bb:bb, bb:sc, sc:bb, sc:sc
residue_residue_contact_map = residue_residue_contact_map(u, chain1, chain2, n_res_1, n_res_2, contact_matrix_bb_bb, contact_matrix_bb_sc, contact_matrix_sc_bb, contact_matrix_sc_sc)
residue_residue_contact_map_bb_bb_normalised = residue_residue_contact_map[0]
residue_residue_contact_map_bb_bb_unnormalised = residue_residue_contact_map[1]
residue_residue_contact_map_bb_sc_normalised = residue_residue_contact_map[2]
residue_residue_contact_map_bb_sc_unnormalised = residue_residue_contact_map[3]
residue_residue_contact_map_sc_bb_normalised = residue_residue_contact_map[4]
residue_residue_contact_map_sc_bb_unnormalised = residue_residue_contact_map[5]
residue_residue_contact_map_sc_sc_normalised = residue_residue_contact_map[6]
residue_residue_contact_map_sc_sc_unnormalised = residue_residue_contact_map[7]

# Save the residue-residue contact map to a .dat file
np.savetxt("residue_residue_contact_map_bb_bb_normalised.dat", residue_residue_contact_map_bb_bb_normalised)
np.savetxt("residue_residue_contact_map_bb_bb_unnormalised.dat", residue_residue_contact_map_bb_bb_unnormalised)
np.savetxt("residue_residue_contact_map_bb_sc_normalised.dat", residue_residue_contact_map_bb_sc_normalised)
np.savetxt("residue_residue_contact_map_bb_sc_unnormalised.dat", residue_residue_contact_map_bb_sc_unnormalised)
np.savetxt("residue_residue_contact_map_sc_bb_normalised.dat", residue_residue_contact_map_sc_bb_normalised)
np.savetxt("residue_residue_contact_map_sc_bb_unnormalised.dat", residue_residue_contact_map_sc_bb_unnormalised)
np.savetxt("residue_residue_contact_map_sc_sc_normalised.dat", residue_residue_contact_map_sc_sc_normalised)
np.savetxt("residue_residue_contact_map_sc_sc_unnormalised.dat", residue_residue_contact_map_sc_sc_unnormalised)

def selected_residue_contact_map(uni, chain1, chain2, n_res_1, n_res_2, contact_matrix_bb_bb, contact_matrix_bb_sc, contact_matrix_sc_bb, contact_matrix_sc_sc):
    """
    Calculate the residue-residue contact map for bb:bb, bb:sc, sc:bb, sc:sc

    Logic:
    ------
    
    
    Returns
    -------
    residue_residue_bb_bb : np.ndarray
        A binary contact map between backbone of chain1 and backbone of chain2
    residue_residue_bb_sc : np.ndarray
        A binary contact map between backbone of chain1 and side chain of chain2
    residue_residue_sc_bb : np.ndarray
        A binary contact map between side chain of chain1 and backbone of chain2
    residue_residue_sc_sc : np.ndarray
        A binary contact map between side chain of chain1 and side chain of chain2
    """
    # Get the sequence of the MUT16 and MUT8 proteins using MDAnalysis
    all_residues = list(uni.residues.resnames)
    MUT16_residues = all_residues[:n_res_1]
    MUT8_residues = all_residues[n_res_1 * chain1 : n_res_1 * chain1 + n_res_2]

    residue_residue_matrix_bb_bb = []
    residue_residue_matrix_bb_sc = []   
    residue_residue_matrix_sc_bb = []
    residue_residue_matrix_sc_sc = []

    # Write a for loop to get the residue-residue contact map for bb:bb, bb:sc, sc:bb, sc:sc

    for i, aa_i in enumerate(MUT16_residues):
        # Initialize the row
        row_bb_bb = []
        row_bb_sc = []
        row_sc_bb = []
        row_sc_sc = []

        for j, aa_j in enumerate(MUT8_residues):
            value_bb_bb = contact_matrix_bb_bb[i, j]
            value_bb_sc = contact_matrix_bb_sc[i, j]
            value_sc_bb = contact_matrix_sc_bb[i, j]
            value_sc_sc = contact_matrix_sc_sc[i, j]

            # Resname and values in one list
            element_bb_bb = [aa_i, aa_j, value_bb_bb]
            element_bb_sc = [aa_i, aa_j, value_bb_sc]
            element_sc_bb = [aa_i, aa_j, value_sc_bb]
            element_sc_sc = [aa_i, aa_j, value_sc_sc]

            # Append the element to the row
            row_bb_bb.append(element_bb_bb)
            row_bb_sc.append(element_bb_sc)
            row_sc_bb.append(element_sc_bb)
            row_sc_sc.append(element_sc_sc)

        # Append the row to the matrix
        residue_residue_matrix_bb_bb.append(row_bb_bb)
        residue_residue_matrix_bb_sc.append(row_bb_sc)
        residue_residue_matrix_sc_bb.append(row_sc_bb)
        residue_residue_matrix_sc_sc.append(row_sc_sc)

    residue_residue_matrix_bb_bb = np.array(residue_residue_matrix_bb_bb)
    residue_residue_matrix_bb_sc = np.array(residue_residue_matrix_bb_sc)
    residue_residue_matrix_sc_bb = np.array(residue_residue_matrix_sc_bb)
    residue_residue_matrix_sc_sc = np.array(residue_residue_matrix_sc_sc)

    # Dictionary for the selected residues only
    aa_dict = {'SER': 0, 'ARG': 1, 'GLY': 2, 'TYR': 3, 'ASP': 4, 'PRO': 5, 'LYS': 6, 'PHE': 7, 'HIS': 8}

    # Helper function to safely get index
    def get_index(residue):
        return aa_dict.get(residue, -1)  # Return -1 if residue is not found

    # Make the residue-residue contact matrix for bb:bb
    aa_matrix_bb_bb = np.zeros((9, 9), dtype=float)
    aa_matrix_count_bb_bb = np.zeros((9, 9), dtype=int)

    for row in residue_residue_matrix_bb_bb:
        for element in row:
            aa_i = element[0]
            aa_j = element[1]
            value = element[2]
            i = get_index(aa_i)
            j = get_index(aa_j)
            if i != -1 and j != -1:  # Only process if both residues are in the dictionary
                aa_matrix_bb_bb[i, j] += float(value)
                aa_matrix_bb_bb[j, i] += float(value)
                aa_matrix_count_bb_bb[i, j] += 1
                aa_matrix_count_bb_bb[j, i] += 1

    aa_matrix_bb_bb_normalised = np.divide(aa_matrix_bb_bb, aa_matrix_count_bb_bb, out=np.zeros_like(aa_matrix_bb_bb), where=aa_matrix_count_bb_bb != 0)
    aa_matrix_bb_bb_unnormalised = aa_matrix_bb_bb

    # Make the residue-residue contact matrix for bb:sc
    aa_matrix_bb_sc = np.zeros((9, 9), dtype=float)
    aa_matrix_count_bb_sc = np.zeros((9, 9), dtype=int)

    for row in residue_residue_matrix_bb_sc:
        for element in row:
            aa_i = element[0]
            aa_j = element[1]
            value = element[2]
            i = get_index(aa_i)
            j = get_index(aa_j)
            if i != -1 and j != -1:
                aa_matrix_bb_sc[i, j] += float(value)
                aa_matrix_bb_sc[j, i] += float(value)
                aa_matrix_count_bb_sc[i, j] += 1
                aa_matrix_count_bb_sc[j, i] += 1

    aa_matrix_bb_sc_normalised = np.divide(aa_matrix_bb_sc, aa_matrix_count_bb_sc, out=np.zeros_like(aa_matrix_bb_sc), where=aa_matrix_count_bb_sc != 0)
    aa_matrix_bb_sc_unnormalised = aa_matrix_bb_sc

    # Make the residue-residue contact matrix for sc:bb
    aa_matrix_sc_bb = np.zeros((9, 9), dtype=float)
    aa_matrix_count_sc_bb = np.zeros((9, 9), dtype=int)

    for row in residue_residue_matrix_sc_bb:
        for element in row:
            aa_i = element[0]
            aa_j = element[1]
            value = element[2]
            i = get_index(aa_i)
            j = get_index(aa_j)
            if i != -1 and j != -1:
                aa_matrix_sc_bb[i, j] += float(value)
                aa_matrix_sc_bb[j, i] += float(value)
                aa_matrix_count_sc_bb[i, j] += 1
                aa_matrix_count_sc_bb[j, i] += 1

    aa_matrix_sc_bb_normalised = np.divide(aa_matrix_sc_bb, aa_matrix_count_sc_bb, out=np.zeros_like(aa_matrix_sc_bb), where=aa_matrix_count_sc_bb != 0)
    aa_matrix_sc_bb_unnormalised = aa_matrix_sc_bb

    # Make the residue-residue contact matrix for sc:sc
    aa_matrix_sc_sc = np.zeros((9, 9), dtype=float)
    aa_matrix_count_sc_sc = np.zeros((9, 9), dtype=int)

    for row in residue_residue_matrix_sc_sc:
        for element in row:
            aa_i = element[0]
            aa_j = element[1]
            value = element[2]
            i = get_index(aa_i)
            j = get_index(aa_j)
            if i != -1 and j != -1:
                aa_matrix_sc_sc[i, j] += float(value)
                aa_matrix_sc_sc[j, i] += float(value)
                aa_matrix_count_sc_sc[i, j] += 1
                aa_matrix_count_sc_sc[j, i] += 1

    aa_matrix_sc_sc_normalised = np.divide(aa_matrix_sc_sc, aa_matrix_count_sc_sc, out=np.zeros_like(aa_matrix_sc_sc), where=aa_matrix_count_sc_sc != 0)
    print (aa_matrix_sc_sc_normalised)
    aa_matrix_sc_sc_unnormalised = aa_matrix_sc_sc

    return (
        aa_matrix_bb_bb_normalised, aa_matrix_bb_bb_unnormalised,
        aa_matrix_bb_sc_normalised, aa_matrix_bb_sc_unnormalised,
        aa_matrix_sc_bb_normalised, aa_matrix_sc_bb_unnormalised,
        aa_matrix_sc_sc_normalised, aa_matrix_sc_sc_unnormalised
    )

# Calculate the residue-residue contact map for bb:bb, bb:sc, sc:bb, sc:sc
residue_residue_contact_map = selected_residue_contact_map(u, chain1, chain2, n_res_1, n_res_2, contact_matrix_bb_bb, contact_matrix_bb_sc, contact_matrix_sc_bb, contact_matrix_sc_sc)

# Unpack the results
residue_residue_contact_map_bb_bb_normalised = residue_residue_contact_map[0]
residue_residue_contact_map_bb_bb_unnormalised = residue_residue_contact_map[1]
residue_residue_contact_map_bb_sc_normalised = residue_residue_contact_map[2]
residue_residue_contact_map_bb_sc_unnormalised = residue_residue_contact_map[3]
residue_residue_contact_map_sc_bb_normalised = residue_residue_contact_map[4]
residue_residue_contact_map_sc_bb_unnormalised = residue_residue_contact_map[5]
residue_residue_contact_map_sc_sc_normalised = residue_residue_contact_map[6]
residue_residue_contact_map_sc_sc_unnormalised = residue_residue_contact_map[7]

# Process the matrix before plotting

#Step1: Transpose the matrix
transposed_residue_residue_contact_map_bb_bb_normalised = np.transpose(residue_residue_contact_map_bb_bb_normalised)
transposed_residue_residue_contact_map_bb_bb_unnormalised = np.transpose(residue_residue_contact_map_bb_bb_unnormalised)
transposed_residue_residue_contact_map_bb_sc_normalised = np.transpose(residue_residue_contact_map_bb_sc_normalised)
transposed_residue_residue_contact_map_bb_sc_unnormalised = np.transpose(residue_residue_contact_map_bb_sc_unnormalised)
transposed_residue_residue_contact_map_sc_bb_normalised = np.transpose(residue_residue_contact_map_sc_bb_normalised)
transposed_residue_residue_contact_map_sc_bb_unnormalised = np.transpose(residue_residue_contact_map_sc_bb_unnormalised)
transposed_residue_residue_contact_map_sc_sc_normalised = np.transpose(residue_residue_contact_map_sc_sc_normalised)
transposed_residue_residue_contact_map_sc_sc_unnormalised = np.transpose(residue_residue_contact_map_sc_sc_unnormalised)

#Step2: Flip the matrix along the y-axis
flipped_residue_residue_contact_map_bb_bb_normalised = np.flipud(transposed_residue_residue_contact_map_bb_bb_normalised)
flipped_residue_residue_contact_map_bb_bb_unnormalised = np.flipud(transposed_residue_residue_contact_map_bb_bb_unnormalised)
flipped_residue_residue_contact_map_bb_sc_normalised = np.flipud(transposed_residue_residue_contact_map_bb_sc_normalised)
flipped_residue_residue_contact_map_bb_sc_unnormalised = np.flipud(transposed_residue_residue_contact_map_bb_sc_unnormalised)
flipped_residue_residue_contact_map_sc_bb_normalised = np.flipud(transposed_residue_residue_contact_map_sc_bb_normalised)
flipped_residue_residue_contact_map_sc_bb_unnormalised = np.flipud(transposed_residue_residue_contact_map_sc_bb_unnormalised)
flipped_residue_residue_contact_map_sc_sc_normalised = np.flipud(transposed_residue_residue_contact_map_sc_sc_normalised)
flipped_residue_residue_contact_map_sc_sc_unnormalised = np.flipud(transposed_residue_residue_contact_map_sc_sc_unnormalised)

#Step3: resultant matrix
resultant_residue_residue_contact_map_bb_bb_normalised = flipped_residue_residue_contact_map_bb_bb_normalised
resultant_residue_residue_contact_map_bb_bb_unnormalised = flipped_residue_residue_contact_map_bb_bb_unnormalised
resultant_residue_residue_contact_map_bb_sc_normalised = flipped_residue_residue_contact_map_bb_sc_normalised
resultant_residue_residue_contact_map_bb_sc_unnormalised = flipped_residue_residue_contact_map_bb_sc_unnormalised
resultant_residue_residue_contact_map_sc_bb_normalised = flipped_residue_residue_contact_map_sc_bb_normalised
resultant_residue_residue_contact_map_sc_bb_unnormalised = flipped_residue_residue_contact_map_sc_bb_unnormalised
resultant_residue_residue_contact_map_sc_sc_normalised = flipped_residue_residue_contact_map_sc_sc_normalised
resultant_residue_residue_contact_map_sc_sc_unnormalised = flipped_residue_residue_contact_map_sc_sc_unnormalised


# Convert the residue-residue contact map to a histogram plot in which the values of bb:bb, bb:sc, sc:bb, sc:sc are plotted on top of each other

# for loop to loop over the selected residues from resultant_residue_residue_contact_map_bb_bb_normalised in x-axis

bar_plot_bb_bb_data = []
bar_plot_bb_sc_data = []
bar_plot_sc_bb_data = []
bar_plot_sc_sc_data = []

for i in range(9):
    for j in range(9-i):
        bar_plot_bb_bb_data.append(resultant_residue_residue_contact_map_bb_bb_normalised[i, j])
        bar_plot_bb_sc_data.append(resultant_residue_residue_contact_map_bb_sc_normalised[i, j])
        bar_plot_sc_bb_data.append(resultant_residue_residue_contact_map_sc_bb_normalised[i, j])
        bar_plot_sc_sc_data.append(resultant_residue_residue_contact_map_sc_sc_normalised[i, j])

# Plot the histogram
fig, ax = plt.subplots()
barWidth = 0.2
r1 = np.arange(9)
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]
r4 = [x + barWidth for x in r3]

plt.bar(r1, bar_plot_bb_bb_data, color='b', width=barWidth, edgecolor='grey', label='bb:bb')
plt.bar(r2, bar_plot_bb_sc_data, color='r', width=barWidth, edgecolor='grey', label='bb:sc')
plt.bar(r3, bar_plot_sc_bb_data, color='g', width=barWidth, edgecolor='grey', label='sc:bb')
plt.bar(r4, bar_plot_sc_sc_data, color='y', width=barWidth, edgecolor='grey', label='sc:sc')
