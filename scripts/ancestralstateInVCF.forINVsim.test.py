#########################################################################################################################################
#																	   #
# This script concatenates the outgroup bases to the end of SNPIDs in a PLINK vcf file							   #
#																		#
# Usage:																			  #
#	python ancestralstate_chr_to_vcf.py  _path_outgroup_chromosomes  vcf_file  output_vcf_file  outlist_removedSNPs  > log				#
#																	 #
# Note that the outgroup_chromosomes' format must be same as those in ~/Desktop/African/GenomeData/humanchimp_hg19/dir_chimp_seq/	#
#########################################################################################################################################


import sys, os, re, gzip


if __name__ == '__main__':

	if sys.argv[1].startswith("std"):
		vcf_obj = sys.stdin
	else:
		vcf_obj = gzip.open(sys.argv[1], 'rb')

	if sys.argv[2].startswith("std"):
		out_vcf_obj = sys.stdout
	else:
		exit ("The 2nd argument must be stdout!")


	fout_mapping_hapID_SV = sys.argv[3]
	fout_mapping_hapID_SV_new = sys.argv[3] + ".new"


	out_snp = []
	
	INFO_AA = False

	while True:
		next_line = vcf_obj.readline()
	
		if next_line.startswith('##'):
			if next_line.startswith('##file'):
				out_vcf_obj.write(next_line)
			else:
				if not INFO_AA:
					INFO_AA = True
					out_vcf_obj.write("##INFO=<ID=AA,Number=A,Type=String,Description=\"Ancestral state based on chimpanzee seq\">\n")
					out_vcf_obj.write("##INFO=<ID=VT,Number=A,Type=String,Description=\"Variant type\">\n")
				else:
					out_vcf_obj.write(next_line)
			continue

		if next_line.startswith("#CHROM"):
			out_vcf_obj.write(next_line)
			l_mapping_hapID_SV = []
			l_mapping_hapID_SV_new = []
			line = next_line.strip().split()
			for sample in line[9:]:
				l_sample = sample.split("_")
				for i, hap in enumerate(l_sample):
					l_mapping_hapID_SV_new.append((sample + "_%s" % (i+1), hap[:2], sample + "_%s" % (i+1)))

					if hap[0] == "D":
						l_mapping_hapID_SV.append((sample + "_%s" % (i+1), "FALSE", sample + "_%s" % (i+1)))
					elif hap[0] == "I":
						l_mapping_hapID_SV.append((sample + "_%s" % (i+1), "TRUE", sample + "_%s" % (i+1)))

			with open(fout_mapping_hapID_SV, "w") as fout:
				fout.write("\t".join(["hapID","SV","orig_hapID"])+"\n")
				for entry in l_mapping_hapID_SV:
					fout.write("\t".join(entry) + "\n")

			with open(fout_mapping_hapID_SV_new, "w") as fout:
				fout.write("\t".join(["hapID","SV","orig_hapID"])+"\n")
				for entry in l_mapping_hapID_SV_new:
					fout.write("\t".join(entry) + "\n")
			continue
			
		if len(next_line) == 0:
			if out_snp != []:
				for entry in out_snp:
					out_vcf_obj.write(entry)
			break

		else:
			snp = next_line.strip().split()
			# skip multiallelic sites
			if len(snp[4]) > 1:
			    continue

			AA_allele = snp[3]

			snp[7] = "AA=" + AA_allele + ";" + "VT=SNP"

			snp = "\t".join(snp)
			out_snp.append(snp+'\n')

			if len(out_snp) >= 1e5:
				for entry in out_snp:
					out_vcf_obj.write(entry)
				out_snp = []
			else:
				continue

