from Bio import SeqIO
import os, sys
import re
import argparse


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='extract FASTA sequences by name')
	parser.add_argument('-i',type=str,dest="infile",required=True,help="Input fasta")
	parser.add_argument('-o',type=str,dest="outfile",required=False, default=None, help="Ouput directory")
	parser.add_argument('-s',type=str,dest="idList",required=False, default=None, help="a list comma-delimited IDs for seq to be removed")
	parser.add_argument('-S',type=str,dest="idFile",required=False, default=None, help="A file with a list of IDs for seq to be removed")
	args = parser.parse_args()

	ll_ids = []

	if args.idList is None and args.idFile is None:
		exit("Error: need -s or -S")
	elif args.idList is not None and args.idFile is not None:
		exit("Error: either -s or -S should be used")
	elif args.idFile is not None:
		with open(args.idFile) as f_ids:
			for ii in f_ids:
				ll_ids.append(ii.strip())
	elif args.idList is not None:
		ll_ids = args.idList.split(",")

	if args.infile in ["-", "std", "stdin"]:
		ifile = sys.stdin
	else:
		ifile = open(args.infile,'rU')

	if args.outfile is not None:
		ofile = open(args.outfile, "w")
	else:
		ofile = sys.stdout

	for record in SeqIO.parse(ifile, "fasta"):
#		if record.name in ll_ids:
		r = [1  for x in ll_ids if re.search(x, record.name)]
		if r != []:
			continue
		else:
			sequence = str(record.seq)
			ofile.write('>'+record.description+'\n'+sequence+'\n')
	
	ofile.close()
