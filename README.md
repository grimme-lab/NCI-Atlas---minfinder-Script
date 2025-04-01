# NCI-Atlas---minfinder-Script
A python3 script to evaluate methods on the NCI Atlas Test sets called "minfinder.py"
The NCI Atlas sets contain: HB300SPXx10, NCIA_D442x10, NCIA_HB375x10, NCIA_IHB100x10, NCIA_SH250x10

For the script tp work, updated .res file versions have to be used in the given benchmark sets
which enable to evaluation of multiple methods in subdirectories, as usually done for benchmark sets.
These .res files are also provided in this repository (make sure to backup the original .res files!)

Descripton of the script and its functionality:

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
