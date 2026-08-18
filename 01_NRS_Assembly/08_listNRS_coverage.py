import numpy

coverage_file = open("SABE1172_UNHESMSV_NRS_dark_freeze_final_coverage.txt", "r")
outfil = open("SABE_1172_UNHESMSV_NRSList_IndividList.txt", "w+")

counter = 0

for line in coverage_file:
 counter += 1
 # Extract list of NRS from header
 if counter == 1:
  samples = line.strip("\n").split("\t")[1:]
  continue
 # Count occurences
 nrs_id = line.split("\t")[0]
 hits = numpy.array([int(h) for h in line.strip("\n").split("\t")[1:]])
 # Get list of RNASeq samples covering the NRS
 presence_index = numpy.nonzero(hits)[0]
 presence_list = [samples[p] for p in presence_index]
 outfil.write(f"{nrs_id}\t{len(presence_list)/len(samples)}\t{','.join(presence_list)}\n")


print("All done, have a nice day!")