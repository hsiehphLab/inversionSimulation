import sys, os, random, math
import msprime
from collections import Counter


#Tsp_p01_p23         /\    
#      1/10 (N_a) I /  \ D  N_a
#                  /    \ 
#Tsp_p0_p1        /\     \
#Tsp_p2_p3       /  \   / \
#               /    \ /   \ 
#              D     I I    D
#             p0    p1 p2    p3
#      1/100(N_a)  1/10(N_a)
#                      1/10(N_a)     N_a
#             |     |  |     |
#             |     _  _     |
#             |       |      |
#             |______P_I_____|
#                     |
#                    P_D


# First we set out the fixed values of the various parameters
chromID = "chr1"
N_a = 6000
mu = 1.25e-8
generation_time = 25

# Times are provided in years, so we convert into generations.
Tsp_p01_p23 = int(sys.argv[1]) / generation_time
Tsp_p0_p1 = int(sys.argv[2]) / generation_time
Tsp_p2_p3 = int(sys.argv[3]) / generation_time


# sample sizes of haplotypes in each pop/deme
#sample_pop0, sample_pop1, sample_pop2, sample_pop3 = [200] * 4
sampleSize = int(sys.argv[4])
inv_freq = float(sys.argv[5])

num_inv = int(sampleSize * inv_freq)
num_inv_sample = round(num_inv/2)
num_direct = sampleSize - num_inv
num_direct_sample = round(num_direct/2)

# Recombination and migration rates
rho = float(sys.argv[6])
m_const = float(sys.argv[7])
seq_length = int(sys.argv[8])

de = msprime.Demography()
de.add_population(name = "P_I", description = "Final INV group", initial_size = 0.1 * N_a)
de.add_population(name = "P_D", description = "Final DIR group", initial_size = N_a)
de.add_population(name = "P0_D", description = "Pop1, DIR", initial_size = 0.01 * N_a)
de.add_population(name = "P1_I", description = "Pop2, INV", initial_size = 0.1 * N_a)
de.add_population(name = "P2_I", description = "Pop3, INV", initial_size = 0.1 * N_a)
de.add_population(name = "P3_D", description = "Pop4, DIR", initial_size = N_a)
de.add_population(name = "Pa_I", description = "Ancestral INV group", initial_size = 0.1 * N_a)
de.add_population(name = "Pa_D", description = "Ancestral DIR group", initial_size = N_a)
de.add_population(name = "P00",  description = "Ancestral group", initial_size = N_a)

de.set_symmetric_migration_rate(["P0_D","P3_D"], m_const)
de.set_symmetric_migration_rate(["P1_I","P2_I"], m_const)

frac_admixI = random.randint(0,10)/10
frac_admixD = random.randint(0,10)/10
de.add_admixture(time=0.00001, derived="P_I", ancestral=["P1_I","P2_I"], proportions = [frac_admixI, 1-frac_admixI])
de.add_admixture(time=0.00001, derived="P_D", ancestral=["P0_D","P3_D"], proportions = [frac_admixD, 1-frac_admixD])
de.add_population_split(time=Tsp_p2_p3, derived=["P2_I","P3_D"], ancestral="Pa_D")
de.add_population_split(time=Tsp_p0_p1, derived=["P0_D","P1_I"], ancestral="Pa_I")
de.add_population_split(time=Tsp_p01_p23, derived=["Pa_I","Pa_D"], ancestral="P00")


l_sampleSize = [ num_inv_sample, num_direct_sample]

sampleIDs = []
for i, v in enumerate(l_sampleSize):
    for ii in range(v):
        if i not in [0]:
            sampleIDs.append("D%s%s_D%s%s" % (i, ii, i, ii))
        else:
            sampleIDs.append("I%s%s_I%s%s" % (i, ii, i, ii))


ts = msprime.sim_ancestry(
        samples=[msprime.SampleSet(round(num_inv_sample), population="P_I", ploidy=2), 
                msprime.SampleSet(round(num_direct_sample), population="P_D", ploidy=2)], demography=de, 
                sequence_length=seq_length, recombination_rate=rho)
mts = msprime.sim_mutations(ts, rate=mu)

with sys.stdout as vcffile:
    mts.write_vcf(vcffile, contig_id = chromID, individual_names=sampleIDs)


