"""
Run to delete all files resulting from compilation and running the srfm.
Deletes:
- rfm outputs
- .o files
- .mod files
- rfm executable
- disort output files
- compiled modules (.so files)
- cache files
- rfm .log file
"""

import os
import shutil

# function to list files
def list_files_walk(start_path='.'):
    all_fnames = []
    for root, dirs, files in os.walk(start_path):
        for file in files:
            all_fnames.append(os.path.join(root, file))
    return all_fnames

# list files
all_fnames = list_files_walk()

#remove files
for file in all_fnames:
    if file.endswith(".asc"):
        os.remove(file)
    elif file.endswith(".o"):
        os.remove(file)
    elif file.endswith(".so"):
        os.remove(file)
    elif file.endswith(".mod"):
        os.remove(file)
    elif file.endswith("rfm"):
        os.remove(file)
    elif file.endswith("code.f"):
        os.remove(file)
    elif file.endswith("INTENSITY.dat"):
        os.remove(file)
    elif file.endswith(".pyc"):
        os.remove(file)
    elif file.endswith("rfm.log"):
        os.remove(file)
    else:
        pass

#directories to remove
dirs = ["./srfm/__pycache__/", "./.idea"]

for d in dirs:
    if os.path.isdir(d):
        shutil.rmtree(d)
    

print("Successfully cleaned srfm.")
