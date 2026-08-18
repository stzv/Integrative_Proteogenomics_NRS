from glob import glob

outfile = open("00a_inputfilenames.txt", "w+")
outfile2 = open("00a_samplenames.txt", "w+")

for fil in glob("00_mzML/*"):
    filename = fil.split("/")[-1]
    samplename = fil.split("/")[-1].replace(".mzML", "")
    outfile.write(f"{{'fileName': '{filename}', 'sampleName': '{samplename}', 'fraction': 1, 'centroiding': True}},\n")
    outfile2.write(f"{{'sampleName': '{samplename}', 'fastaName': '{samplename}', 'paramName': 'normalSearchMSFragger'}},\n")