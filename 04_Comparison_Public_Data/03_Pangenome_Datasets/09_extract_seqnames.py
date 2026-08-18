from glob import glob
import gzip
import sys

setname = sys.argv[1]

files = glob(f"04_mapped_SABE1172_UNHESMSV_{setname}.sam.gz") + glob(f"04_partial_SABE1172_UNHESMSV_{setname}.sam.gz")

list_output = open(f"09_SABE1172_UNHESMSV_{setname}_aligned_list.txt", "w")
#list_output.write("ctg,length_bp,chromosome,position\n")
entries = []

for fil in files:
    print(f"Processing {fil}")
    f = gzip.open(fil,'rb')
    f = f.readlines()
    for line in f:
        string = line.split(b"\t")
        str_length = string[0].split(b"_")
        length = (int(str_length[3]) - int(str_length[2]) + 1)
        chrom = string[2].decode()
        pos = string[3].decode()
        #print(string[0].decode("utf-8"), length)
        entries.append(f'{string[0].decode("utf-8")},{length},{chrom},{pos}')

if len(entries) == len(set(entries)):
    print("There are no duplicate IDs")
else:
    print("There are duplicate IDs, cancelling the run")


for item in set(entries):
    list_output.write(f"{item}\n")


list_output.close()


    
