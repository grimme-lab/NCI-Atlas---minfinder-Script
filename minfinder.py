# This script was created by Benedikt Baedorf to evaluate different Rezak benchmark sets.
# it should work for: HB300SPXx10, NCIA_D442x10, NCIA_HB375x10, NCIA_IHB100x10 and NCIA_SH250x10

import os
import re
import sys
import argparse
import numpy as np
from scipy.interpolate import CubicSpline

#----------------------------- HELP PRINTOUT --------------------------------------------------------------------------------

def print_help_page():
    print("""
    Usage: python3 minfinder.py [method_name]
        or python3 minfinder.py ref 
        or python3 minfinder.py ref [method_name]
        or python3 minfinder.py [method_name] -nosub 

    Description:
      This script evaluates the ".res" file and evaluates the given PES scan data (either of ther reference or of a given method to evaluate).
      By fitting a spline function to the given datapoints, the script determines the interpolated minimum of the PES.
      It then computes the scaling factor for that given minimum. After determining the closest intermolecular distances in the dimer structures 
      (for scaling factor 1.0), the script then applies the determined minimum scaling factor to obtain a interpolated distance at which
      the PES minimum for a given method would lie [unit: angstrom].

      If 'ref' is provided as the method, the script extracts reference minimum values from the ".res" file.
      Otherwise, it processes the method-specific directories that should lie as subdirectories in each system folder.
      If the "-nosub" flag is given, the script does not look for subdirectories but just the energy files in the system folders.
      The script can ready 'energy' files and 'orca.out' files.
      If you want to reduce the amount of systems to evalaute, make a reduced version of the ".res" file. (Make sure to backup your .res file)

    Options:
      method_name       Evaluate the given method (e.g. g-xtb, gfn2-xtb, ...)   --> Prinout: systemid intermol_dist_PES_min_method    [in angstrom]
      ref               Evaluate the minima of the reference method             --> Prinout: systemid intermol_dist_PES_min_reference [in angstrom]
      ref method_name   Returns the shortest intermolecular distance for        --> Prinout: intermol_dist_PES_min_reference intermol_dist_PES_min_method [in angstrom]
                        each system for the reference and the method
      -nosub            If the method results are not in subdirectories        
      -help             Show this help message and exit

    """)
    sys.exit(0)

if "-help" in sys.argv or "--help" in sys.argv or "-h" in sys.argv:
    print_help_page()

#---------------------------- EXTRACTING systemid, scaling factor and folder  ------------------------------------------------

def find_reference_minimum(res_file=".res", method_name=None, use_subfolder=True):
    if not os.path.exists(res_file):
        print(f"Error: '{res_file}' not found.")
        return {}, {}

    system_id_pattern = re.compile(r"([\d\.]+)_")  # System ID pattern
    scaling_factor_pattern = re.compile(r"_(\d[\d\.]{2,3})")  # Scaling factor pattern
    folder_pattern = re.compile(r"\{([^}]+)\}")  # Extract folder names in brackets
    reference_data = {}  # Dictionary to store extracted data
    folder_data = {}     # Stores folders for scaling factor 1.0

    with open(res_file, "r") as f:
        for line in f:
            system_id_match = system_id_pattern.search(line)
            if system_id_match:
                system_id = system_id_match.group(1)
                scaling_factor_match = scaling_factor_pattern.search(line)
                
                if scaling_factor_match:
                    raw_scaling_factor = scaling_factor_match.group(1)
                    scaling_factor = normalize_scaling_factor(raw_scaling_factor)
                    
                    if scaling_factor is not None:
                        folder_match = folder_pattern.search(line)
                        if folder_match:
                            folder_str = folder_match.group(1)
                            folders = folder_str.split(",")
                            
                            for folder in folders:
                                folder = folder.strip()
                                
                                # Check if the system_id and scaling_factor match the folder name
                                if system_id in folder and (f"_{str(scaling_factor)}" in folder or f"_{str(raw_scaling_factor)}" in folder):
                                    
                                    # Store folder data for scaling factor 1.0
                                    if scaling_factor == 1.0 and len(folders) == 3:
                                        folder_data[system_id] = tuple(folders)

                                    if method_name == "ref":  # Extract energy from the last value in the line
                                        try:
                                            energy = float(line.split()[-1])  # Extract last entry as energy
                                            if system_id not in reference_data:
                                                reference_data[system_id] = []
                                            reference_data[system_id].append((scaling_factor, energy))
                                        except ValueError:
                                            print(f"Warning: Could not extract energy for system {system_id} in .res file")
                                    else:
                                        # Adjust the method folder logic based on `use_subfolder`
                                        if use_subfolder:
                                            method_folder = os.path.join(str(folder), str(method_name))
                                            #print(method_folder)
                                        else:
                                            method_folder = os.path.join(str(folder))  # Directly use the folder itself if no subfolder

                                        #print(f"Looking in folder: {method_folder}")
                                        
                                        # Check if the directory exists
                                        if os.path.isdir(method_folder):
                                            energy = extract_energy(method_folder, method_name, use_subfolder)
                                            if energy is not None:
                                                if system_id not in reference_data:
                                                    reference_data[system_id] = []
                                                reference_data[system_id].append((scaling_factor, energy))
                                        else:
                                            print("The folder", method_folder,"does not exist. Exiting ...")
                                            sys.exit()

    return reference_data, folder_data

#-------------------------------------- EXTRACT AND ADJUST SCALING FACTOR ---------------------------------------------

def extract_scaling_factor(name):
    match = re.search(r'_(\d)([\d\.]{2,3})', name)
    if match:
        raw_value = match.group(1) + match.group(2)
        return normalize_scaling_factor(raw_value)
    return None

def normalize_scaling_factor(raw_value):
    """
    Convert raw scaling factor strings like '080' into proper float values (0.80).
    Ensures a correct decimal placement.
    """
    if isinstance(raw_value, str):
        raw_value = raw_value.replace("_", ".")  # Replace underscores with dots if needed
        try:
            if "." not in raw_value and len(raw_value) > 1:
                raw_value = raw_value[:1] + "." + raw_value[1:]
            return float(raw_value)
        except ValueError:
            print(f"Warning: Could not normalize scaling factor '{raw_value}'")
            return None
    return None

#--------------------------------------- EXTRACT ENERGY --------- -----------------------------------------------

def extract_energy(folder, method_name, use_subfolder=True):
    """
    Extract energy from the given folder.

    Parameters:
    - folder (str): The directory to look in.
    - method_name (str): The method (e.g., "g-xtb").
    - use_subfolder (bool): Whether to look inside method subfolders or directly in the folder.

    Returns:
    - float or None: Extracted energy value.
    """
    energy_file = os.path.join(folder, "energy")
    orca_out = os.path.join(folder, "orca.out")

    if os.path.isfile(orca_out):
        with open(orca_out, "r") as f:
            for line in f:
                if "FINAL SINGLE POINT ENERGY" in line:
                    try:
                        return float(line.split()[-1])
                    except (IndexError, ValueError):
                        print(f"Warning: Could not extract energy from {orca_out}")
                        return None
    elif os.path.isfile(energy_file):
        with open(energy_file, "r") as f:
            lines = f.readlines()
            if len(lines) > 1:
                try:
                    return float(lines[1].split()[1])
                except (IndexError, ValueError):
                    print(f"Warning: Invalid energy value in {energy_file}")
                    return None
    return None

#--------------------------------------- FIND SCALING FACTOR FOR SPLINE MINIMUM -----------------------------------

def find_spline_minimum(energy_data):
    """
    Finds the scaling factor (x-value) at which the spline function has its minimum
    for a given system_id dictionary.

    Parameters:
    - energy_data (dict): Dictionary with a single system_id as the key and a list of
      (scaling_factor, energy) tuples as the value.

    Returns:
    - float: The x-value corresponding to the minimum of the spline for the given system_id.
    """
    if not energy_data:
        return None  # No data to process

    # Extract system_id and data
    system_id, data = list(energy_data.items())[0]

    if not data:
        return None  # No valid data

    data.sort()  # Ensure data is sorted by scaling factor
    x = np.array([point[0] for point in data])  # Scaling factors
    y = np.array([point[1] for point in data])  # Energies

    if len(x) != len(set(x)):  # Check for duplicates
        print(f"Error: Duplicate scaling factor found for system {system_id}")
        return None

    spline = CubicSpline(x, y)

    # Find the minimum by evaluating the derivative
    derivative = spline.derivative()
    critical_points = derivative.roots()

    # Filter critical points within the data range
    valid_critical_points = [pt for pt in critical_points if x[0] <= pt <= x[-1]]

    if valid_critical_points:
        min_x = round(min(valid_critical_points, key=lambda pt: spline(pt)), 3)
        return min_x  # Return only the min_x value
    else:
        print(f"{system_id} No valid minimum found")
        return None

#--------------------------------------- FIND SHORTEST DISTANCE FOR SCALE FACTOR 1.0 -------------------------------

def find_min_distance_in_dimer(monomerA_file, monomerB_file, dimer_file):
    """Finds the minimum distance between monomers in a dimer structure and returns it in Angstroms."""
    def read_coord_file(filename):
        """Reads a coord file and extracts atomic coordinates."""
        with open(filename, 'r') as f:
            lines = f.readlines()

        coords = []
        read = False
        for line in lines:
            line = line.strip()
            if line == "$coord":
                read = True
                continue
            elif line == "$end":
                break
            elif read and line:
                parts = line.split()
                if len(parts) == 4:  # x, y, z, atom_type
                    x, y, z = map(float, parts[:3])
                    coords.append((x, y, z))

        return np.array(coords)

    def calculate_min_distance(coords1, coords2):
        """Calculates the minimum distance between two sets of coordinates."""
        min_distance = float('inf')
        for atom1 in coords1:
            for atom2 in coords2:
                distance = np.linalg.norm(np.array(atom1) - np.array(atom2))
                if distance < min_distance:
                    min_distance = distance
        return min_distance * 0.529177  # Convert from Bohr to Angstrom

    coordsA = read_coord_file(monomerA_file)
    coordsB = read_coord_file(monomerB_file)
    coordsDimer = read_coord_file(dimer_file)

    # Extract monomer coordinates from dimer file based on assumed ordering
    coordsDimerA = coordsDimer[:len(coordsA)]  # First part is monomer A
    coordsDimerB = coordsDimer[len(coordsA):]  # Second part is monomer B

    return calculate_min_distance(coordsDimerA, coordsDimerB)

#====================================================== MAIN ============================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract and analyze energy values for each system.")
    parser.add_argument("method_name", nargs="+", type=str, help="Method name(s) (e.g., ref, g-xtb, etc.)")
    parser.add_argument("-nosub", action="store_true", help="Disable subfolder usage (use energy directly in folder).")
    args = parser.parse_args()

    # A dictionary to hold results for each method independently
    method_results = {}
    method_scaled_distances = {}  # Dictionary to store scaled distances for each method

    # Determine whether to use subfolders or not
    use_subfolder = not args.nosub  # If -nosub is provided, we disable subfolder use

    # Process each method name passed in args.method_name
    for method in args.method_name:
        # Initialize separate dictionaries for each method
        method_results[method] = {}
        method_scaled_distances[method] = {}

        if method == "ref":
            # Handle the "ref" method separately
            reference_data, folder_data = find_reference_minimum(res_file=".res", method_name="ref", use_subfolder=use_subfolder)

            # Now correctly populate the method_results with reference data
            for system_id, data in reference_data.items():  # 'data' contains list of (scaling_factor, energy) tuples
                for scaling_factor, energy in data:
                    if energy is not None:
                        if system_id not in method_results[method]:
                            method_results[method][system_id] = []  # Initialize the list if not already
                        method_results[method][system_id].append((scaling_factor, energy))
        else:
            # Handle other methods, including "g-xtb"
            reference_data, folder_data = find_reference_minimum(res_file=".res", method_name=method, use_subfolder=use_subfolder)

            # Now correctly populate the method_results with extracted data
            for system_id, data in reference_data.items():  # 'data' contains (scaling_factor, energy) tuples
                for scaling_factor, energy in data:
                    if energy is not None:
                        if system_id not in method_results[method]:
                            method_results[method][system_id] = []  # Initialize if not already
                        method_results[method][system_id].append((scaling_factor, energy))  # Same as "ref"

    # Calculate minimum distances
    distance_data = {}  # Dictionary to store system_id and corresponding min distance
    for system_id, folders in folder_data.items():
        if len(folders) == 3:  # Ensure there are exactly 3 folders
            monomerA_path = os.path.join(folders[1], "coord")  # Monomer A (index 1)
            monomerB_path = os.path.join(folders[2], "coord")  # Monomer B (index 2)
            dimer_path = os.path.join(folders[0], "coord")     # Dimer (index 0)

            min_distance = find_min_distance_in_dimer(monomerA_path, monomerB_path, dimer_path)
            distance_data[system_id] = min_distance  # Store result

    # Compute scaled minimum distances for each method and system
    for method in args.method_name:
        for system_id, data in method_results[method].items():  # Iterate through the dictionary for each method
            min_x = find_spline_minimum({system_id: data})  # Find the scaling factor for the minimum

            if min_x is not None and system_id in distance_data:  # Ensure valid scaling factor and distance data
                scaled_distance = min_x * distance_data[system_id]  # Compute scaled distance
                method_scaled_distances[method][system_id] = scaled_distance  # Store result

#-------------------------------------- PRINT OUTS -----------------------------------------------------------------
# Output the results based on the number of methods provided
if len(args.method_name) == 1:
    # For a single method, print system_id and the scaled distance for that method
    for system_id in method_scaled_distances[args.method_name[0]].keys():
        scaled_dist = method_scaled_distances[args.method_name[0]].get(system_id, "N/A")  # Get distance or "N/A"
        print(system_id, f"{scaled_dist:.3f}")  # Print system_id and scaled distance with 3 decimal places

else:
    # For multiple methods, print only the scaled distances, without system_id
    system_ids = list(method_scaled_distances[args.method_name[0]].keys())  # Get list of system_ids

    for system_id in system_ids:
        for method in args.method_name:
            scaled_dist = method_scaled_distances[method].get(system_id, "N/A")  # Get distance or "N/A"
            print(f"{scaled_dist:.3f}", end=" ")  # Print the scaled distance for the current method
        print()  # Move to the next line after printing distances for all methods

