import subprocess
import os

for line in open("PXD020192_README.txt", "r"):
    ID, NAME, URI, TYPE, MAPPINGS = line.rstrip("\n").split("\t")
    # Skip if file already donwloaded
    if os.path.isfile(NAME) == True:
        continue
    # Download the file
    if URI.endswith(".raw"):
        subprocess.run(f"wget {URI}", shell = True)

for line in open("PXD026370 _checksum.txt", "r"):
    if line.startswith("#"):
        continue
    #
    fil, md5 = line.strip("\n").split("\t")
    fil_id = fil.split("\\")[-1]
    URI = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/08/PXD026370/" + fil_id
    # Skip if file already downloaded
    if os.path.isfile(fil_id) == True:
        continue
    if URI.endswith(".raw"):
        subprocess.run(f"wget {URI}", shell = True)