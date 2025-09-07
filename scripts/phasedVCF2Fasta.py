
## This script is used to calculate the observed divergence between human and chimp using a sample of human genomes. The divergence is taken by averaging the divergence estimates between the chimp ref and each of the samples (haploid sequences) according to the cPickle and TPED file.
##
## Usage:
##		python new_calculate_divergence_humanchimp_seq_v3_average_div_between_chimp_and_TPED.py  _path_chimpRef  _path_humanRef   _filename_cPickle_biteMap  _option_calculation  __fname_TPEDfile  >  chimp_human_divergence.dat
##
## The input file '_filename_cPickle_biteMap' is typically based on called positions across all samples.
## This is a correct calculation for sequence divergence based on Nei 1987 (equation 5.3 and 10.20).

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import random, os, re, sys, time, gzip, vcf, pysam, copy
import numpy as np
import argparse
import pickle as cPickle
from pyfaidx import Fasta

parser = argparse.ArgumentParser()
parser.add_argument("--vcf")
parser.add_argument("--fout")
parser.add_argument("--locus", nargs="?", const=1, default=None)
parser.add_argument("--addAAseq", nargs="?", const=1, default=None)
parser.add_argument("--reference", nargs="?", const=1, default=None)
args = parser.parse_args()


def func_overlap(x, y):
	return range(max(x[0], y[0]), min(x[-1], y[-1])+1)

reference_seq = Fasta(args.reference)
chromID, start_end = args.locus.split(":")
pos_start, pos_end = [int(x) for x in start_end.split("-")]

vcf_reader = vcf.Reader(filename=args.vcf)
vcf_records = vcf_reader.fetch(chromID, pos_start, pos_end)
total_number_haploidSamples = 2 * len(vcf_records.samples)

list_sampleID_seqs = []

if args.addAAseq:
	AA_seq = list(reference_seq[chromID][ (pos_start-1) : pos_end ].seq)


for sample in vcf_records.samples:
	human_seq1 = list(reference_seq[chromID][ (pos_start-1) : pos_end ].seq)
	human_seq2 = list(reference_seq[chromID][ (pos_start-1) : pos_end ].seq)
	list_sampleID_seqs.append([sample, human_seq1, human_seq2])

list_outVCF = []
for record in vcf_records:

	lookup_REF_ALT = [str(record.REF[0])]
	lookup_REF_ALT.extend([str(x) for x in record.ALT])
	if args.addAAseq:
		AA_allele = record.INFO["AA"][0]
		try:
			AA_seq[record.POS - pos_start] = AA_allele

		except IndexError:
			print("AAseq")
			print(AA_allele, record.POS, pos_start)
			print(record.POS - pos_start)
			print(len(AA_seq))
			exit("IndexError AAseq")

		
	for idx, sample in enumerate(record.samples):
		if sample["GT"] != ".":
			try:
			    alleles = [lookup_REF_ALT[int(x)] if x != "." else "N" for x in re.split("\/|\|", sample["GT"])]
			except IndexError:
			    print (record.POS)
			    exit("WTF")
		else:
			alleles = ["N", "N"]

		try:
			list_sampleID_seqs[idx][1][record.POS - pos_start] = alleles[0]
			list_sampleID_seqs[idx][2][record.POS - pos_start] = alleles[1]

		except IndexError:
			print(sample)
			print(lookup_REF_ALT)
			print(alleles, record.POS, pos_start)
			print(record.POS - pos_start)
			print(len(list_sampleID_seqs[idx][1]))
			print(list_sampleID_seqs[idx][1][record.POS - pos_start])
			print(len(list_sampleID_seqs[idx][1]))
			exit("IndexError")


out_list = []
if args.addAAseq:
	out_list.append(">%s\n" % args.addAAseq)
	seq1 = ''.join(AA_seq) + "\n"
	out_list.append(seq1)
	
for sample in list_sampleID_seqs:
	sampleID = sample[0]

	for seq in range(1,3):
		hapID = sampleID + "_" + str(seq)
		out_list.append(">" + hapID + "\n")
		seq1 = ''.join(sample[seq]) + "\n"
		out_list.append(seq1)


if args.fout in ["-", "std", "stdout"]:
	fout = sys.stdout
else:
	fout = open(args.fout, "w")

for entry in out_list:
	fout.write(entry)

