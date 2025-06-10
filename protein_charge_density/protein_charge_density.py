import MDAnalysis as mda
import numpy as np
import os

def calculate_protein_charge_count(u, z_slice=1):
    """
    Calculate the number of protein cations (ARG, LYS) and anions (ASP, GLU) in each z-slice.
    """
    box_z = u.dimensions[2]

    # Initialize arrays to store the number of cationic and anionic residues in each slice
    protein_cation = np.zeros(int(box_z // z_slice))  # Positive protein charges (ARG, LYS)
    protein_anion = np.zeros(int(box_z // z_slice))   # Negative protein charges (ASP, GLU)

    # Precompute atom selections for cationic and anionic residues
    arg_lys_atoms = u.select_atoms('resname ARG LYS')  # Positively charged residues
    asp_glu_atoms = u.select_atoms('resname ASP GLU')  # Negatively charged residues

    n_frames = 0
    for ts in u.trajectory:
        n_frames += 1
        # Loop over the z-axis in slices of size z_slice
        for i in range(0, int(box_z), z_slice):
            z_start, z_end = i, i + z_slice
            slice_index = i // z_slice  # Corresponding index for the slice

            # Count the number of cationic (ARG, LYS) and anionic (ASP, GLU) residues in the current slice
            protein_cation[slice_index] += np.sum((arg_lys_atoms.positions[:, 2] >= z_start) & (arg_lys_atoms.positions[:, 2] < z_end))
            protein_anion[slice_index] += np.sum((asp_glu_atoms.positions[:, 2] >= z_start) & (asp_glu_atoms.positions[:, 2] < z_end))

    # Normalize by the number of frames to get the average number of atoms in each slice
    protein_cation /= n_frames
    protein_anion /= n_frames

    return protein_cation, protein_anion

def main(filename):
    topology = "/mnt/home/gkumar/Analysis_new2/frame100.gro"

    if os.path.exists(filename) and os.path.isfile(filename):
        u = mda.Universe(topology, filename)
        print(u, flush=True)

        # Calculate the average number of cations and anions in each z-slice
        protein_cation_count, protein_anion_count = calculate_protein_charge_count(u, z_slice=1)

        # Save the atom counts to a file
        np.savetxt('protein_charge_count.txt',
                   np.column_stack((protein_cation_count, protein_anion_count)),
                   header='Protein_Cation_Count (ARG,LYS) Protein_Anion_Count (ASP,GLU)',
                   fmt='%10.5f')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Calculate the number of protein cations and anions in each z-slice.')
    parser.add_argument('filename', type=str, help='The path to the trajectory file (e.g., traj.xtc, traj.dcd).')
    args = parser.parse_args()

    main(args.filename)

