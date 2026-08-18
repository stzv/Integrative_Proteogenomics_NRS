from glob import glob

all_files = glob("brazil/20210817_STAR/results/*/*.bam")

domain_dict = dict()

for fil in all_files:
 domain = fil.split("/")[-2]
 rna = fil.split("/")[-1].strip(".bam")
 #
 if domain not in domain_dict.keys():
  domain_dict[domain] = set()
 #
 domain_dict[domain].add(rna)

outfil = open("03_SABE_1172_UNHESMSV_ACCESSIDs_RNAIDs.txt", "w")

for domain, rnaids in domain_dict.items():
 outfil.write(f"{domain},{len(rnaids)},{','.join(rnaids)}\n")