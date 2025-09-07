import sys, os, re, gzip



def vcf_chrom_iter(vcf_filename, RECODE, RECODE_AA):

	if vcf_filename in ["STDIN","stdin","-"]:
		vcfObj = sys.stdin
	else:
		vcfObj = open(vcf_filename, 'r')

	chrom = {}
	curr_chrom = 0

	while True:
		line = vcfObj.readline()
		if line == "":
			break
		if line.startswith('#'):
			if line.startswith('#CHROM'):
				vcfHeader = line.strip().split()
				yield (None, vcfHeader)
			continue

		l_locus = line.strip().split()
		tmp_snp = [l_locus[0], l_locus[2], "0", l_locus[1]]
		REF = l_locus[3]
		ALT = l_locus[4].split(",")
		alleles = [l_locus[3]]
		alleles.extend(ALT)
		pos = int(l_locus[1])

		genos = l_locus[9:]
		for g in genos:
			g = g.split(":")[0]
			if RECODE:
				try:
					tmp_snp.extend([alleles[int(x)] for x in re.split("[|/]", g)])
				except ValueError:
					exit (g)
			else:
				tmp_genos = [x for x in re.split("[|/]", g)]
				if RECODE_AA:
					AA = l_locus[7].split(";")[0].split("=")[1]
					AA_idx = alleles.index(AA)
					if AA_idx == 1:
						tmp_genos = [ "1" if x == "0" else "0" if x == "1" else "." for x in tmp_genos]
				tmp_snp.extend(tmp_genos)


		try:
			chromID = int(l_locus[0])
		except ValueError:
			if l_locus[0].startswith("chr"):
				try:
					chromID = int(l_locus[0][3:])
				except ValueError:
					chromID = l_locus[0][3:]
			else:
				chromID = l_locus[0]

		if chromID == curr_chrom:
			chrom[pos] = ' '.join(tmp_snp) + '\n'
		elif chromID != curr_chrom:
			yield  (curr_chrom, chrom)
			chrom = {}
			curr_chrom = chromID
			chrom[pos] = ' '.join(tmp_snp) + '\n'

	yield  (curr_chrom, chrom)
			

if __name__ == "__main__":

	fname_vcf = sys.argv[1]

	outPrefix = sys.argv[2]

	# recode numeric to base genotypes
	try:
		recode = eval(sys.argv[3])
	except IndexError:
		recode = None

	try:
		recodeAA = eval(sys.argv[4])
	except IndexError:
		recodeAA = None


	foutname_tped = outPrefix + '.tped'
	if os.path.isfile(foutname_tped):
		os.remove(foutname_tped)
	fout_tped = open(foutname_tped, 'w')

	foutname_tfam = outPrefix + '.tfam'
	if os.path.isfile(foutname_tfam):
		os.remove(foutname_tfam)
	fout_tfam = open(foutname_tfam, 'w')

	for chromosome, snps_chrom_dict in vcf_chrom_iter(fname_vcf, recode, recodeAA):
		if chromosome == None:
			if snps_chrom_dict[0] == "#CHROM":
				list_tfam = []
				IDs = snps_chrom_dict[9:]
				print (IDs)
				list_tfam = [IDs, IDs, [0]*len(IDs), [0]*len(IDs), [0]*len(IDs), [-9]*len(IDs)]
				for row in zip(*list_tfam):
					fout_tfam.write(' '.join([str(x) for x in row]) + '\n')
				fout_tfam.close()

		else:
			if chromosome == 0:
				continue
			for site in sorted(snps_chrom_dict):
				fout_tped.write(snps_chrom_dict[site])

	fout_tped.close()


