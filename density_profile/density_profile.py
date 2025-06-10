import MDAnalysis as mda
import numpy as np
import os

def calculate_atom_counts(u, z_slice=1):
    # Fixed dimensions of the box
    box_z = u.dimensions[2]

    # Initialize arrays to store the atom counts, now with more slices due to smaller bin size
    atom_counts_mut16 = np.zeros(int(box_z // z_slice))
    atom_counts_mut8 = np.zeros(int(box_z // z_slice))
    atom_counts_water = np.zeros(int(box_z // z_slice))
    atom_counts_cl_ions = np.zeros(int(box_z // z_slice))
    atom_counts_na_ions = np.zeros(int(box_z // z_slice))

    # Precompute atom selections
    mut16_atoms = u.select_atoms('index 0:220399')
    mut8_atoms = u.select_atoms('index 220400:228449')
    water_atoms = u.select_atoms('resname SOL')
    cl_ions = u.select_atoms('resname CL')
    na_ions = u.select_atoms('resname NA')

    # Atom positions stored for each atom group
    atom_data = {
        'mut16': {'count_array': atom_counts_mut16, 'positions': mut16_atoms.positions[:, 2]},
        'mut8': {'count_array': atom_counts_mut8, 'positions': mut8_atoms.positions[:, 2]},
        'water': {'count_array': atom_counts_water, 'positions': water_atoms.positions[:, 2]},
        'cl': {'count_array': atom_counts_cl_ions, 'positions': cl_ions.positions[:, 2]},
        'na': {'count_array': atom_counts_na_ions, 'positions': na_ions.positions[:, 2]},
    }

    n_frames = 0
    for ts in u.trajectory:
        n_frames += 1
        # Loop over the z-axis in slices of size z_slice = 1 Å
        for i in range(0, int(box_z), z_slice):
            z_start, z_end = i, i + z_slice
            slice_index = i // z_slice  # Corresponding index for the slice

            for key, value in atom_data.items():
                # Count the number of atoms in the current z_slice
                num_atoms_in_slice = np.sum((value['positions'] >= z_start) & (value['positions'] < z_end))
                # Add the count to the corresponding slice
                value['count_array'][slice_index] += num_atoms_in_slice

    # Return the average number of atoms per bin, divided by the number of frames
    return [value['count_array'] / n_frames for value in atom_data.values()]

def main(filename):
    topology = "/mnt/home/gkumar/Analysis_new2/frame100.gro"

    if os.path.exists(filename) and os.path.isfile(filename):
        u = mda.Universe(topology, filename)
        print(u, flush=True)

        avg_atom_counts = calculate_atom_counts(u, z_slice=1)
        np.savetxt('atom_count_profile_1.txt',
                   np.column_stack(avg_atom_counts),
                   header='Mut16_atoms  Mut8_atoms  Water_atoms  Cl_atoms  Na_atoms',
                   fmt='%10.5f')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Calculate atom counts in each z-bin from a molecular dynamics simulation.')
    parser.add_argument('filename', type=str, help='The path to the trajectory file (e.g., traj.xtc, traj.dcd).')
    args = parser.parse_args()

    main(args.filename)

