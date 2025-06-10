import os
import numpy as np

def get_density_profile_files(base_dir):
    """
    Collect all the 'density_profile.txt' files from the RUN/CLONE/results directories.
    """
    density_files = []
    
    # Loop through the RUN, CLONE, and result folders
    for run in range(40):
        for clone in range(100):
            for result in range(100):
                density_profile_path = os.path.join(base_dir, f"RUN{run}", f"CLONE{clone}", f"results{result}", 'atom_count_profile_1.txt')
                
                # Check if the density_profile.txt file exists
                if os.path.exists(density_profile_path):
                    density_files.append(density_profile_path)
    
    return density_files

def sum_and_average_density_profiles(base_dir, output_file):
    """
    Sum, average, and calculate the standard error of the mean (SEM) of the first 300 rows of all density_profile.txt files.
    """
    density_files = get_density_profile_files(base_dir)
    
    if not density_files:
        print("No density_profile.txt files found.")
        return
    
    # Initialize sum array, sum of squares array, and counter
    sum_array = None
    sum_of_squares = None
    count = 0
    
    # Loop through the files and sum the first 300 rows of the contents
    for file in density_files:
        data = np.loadtxt(file)
        
        # Trim to the first 300 rows
        if data.shape[0] > 300:
            data = data[:300, :]  # Use only the first 300 rows

        # Initialize sum_array and sum_of_squares on the first file
        if sum_array is None:
            sum_array = np.zeros_like(data)  # Initialize sum_array with the same shape as the data
            sum_of_squares = np.zeros_like(data)  # Initialize sum_of_squares for calculating variance
        
        # Make sure all data arrays have the same number of rows (300)
        if data.shape[0] == 300:
            sum_array += data  # Sum the columns
            sum_of_squares += data**2  # Sum the squares of the columns for variance calculation
            count += 1

    # Calculate average and standard error of the mean
    if count > 0:
        average_array = sum_array / count
        
        # Calculate variance and standard deviation
        variance_array = (sum_of_squares / count) - (average_array**2)
        
        # Ensure variance is non-negative (clip negative values to zero)
        variance_array = np.clip(variance_array, 0, None)
        
        stddev_array = np.sqrt(variance_array)
        
        # Calculate the standard error of the mean (SEM)
        sem_array = stddev_array / np.sqrt(count)
        
        # Save the averaged density profile and SEM
        output_data = np.column_stack((average_array, sem_array))
        np.savetxt(output_file, output_data, fmt='%10.5f', header='Averaged Density Profiles and SEM', comments='')
        print(f"Averaged density profile and SEM saved to: {output_file}")
    else:
        print("No files with 300 rows were found.")

# Define the base directory and output file
base_directory = "/mnt/ceph/users/gkumar/data_new"  # Replace with your actual directory path
output_file = "/mnt/home/gkumar/Mut16_Mut8_atomistic_project_analysis/density_profile_average/averaged_number_density_profile_and_sem.txt"

# Run the function to sum, average, and calculate the SEM of the density profiles
sum_and_average_density_profiles(base_directory, output_file)

