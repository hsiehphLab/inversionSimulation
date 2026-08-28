import sys
import msprime

# Tsp_p0_p1                  /\   
#                           D /  \ I
#                            /    \ 
#                           /      \ 
#                          D        I
#                         p0        p1


# First we set out the fixed values of the various parameters
N_a = 6000
mu = 1.25e-8
generation_time = 25

# Times are provided in years, so we convert into generations.
Tsp_p0_p1 = int(sys.argv[1]) / generation_time

# sample sizes of haplotypes in each pop/deme
sample_pop0, sample_pop1 = eval(sys.argv[2])

diploid_sample_pop0 = sample_pop0 // 2
diploid_sample_pop1 = sample_pop1 // 2

# Recombination and migration rates
rho = float(sys.argv[3])
m_const = float(sys.argv[4]) # Parsed but unused in the original migration matrix
seq_length = int(sys.argv[5])
chromID = sys.argv[6]

# Generate sample names
haploIDs = []
for i, v in enumerate([sample_pop0, sample_pop1]):
    for ii in range(v):
        if i != 1:
            haploIDs.append(f"D{i}{ii}")
        else:
            haploIDs.append(f"I{i}{ii}")

sampleIDs = []
for i in range(0, len(haploIDs), 2):
    sampleIDs.append(f"{haploIDs[i]}_{haploIDs[i+1]}")

# 1. Define the Demography
demography = msprime.Demography()

# Initially, 0=pop0, 1=pop1
demography.add_population(name="pop0", initial_size=N_a)
demography.add_population(name="pop1", initial_size=N_a / 100)

# Merge pop1 into pop0 backwards in time at Tsp_p0_p1
# In forwards-time terminology, pop1 is derived from an ancestral pop0
demography.add_population_split(time=Tsp_p0_p1, derived=["pop1"], ancestral="pop0")


# 2. Simulate Ancestry (Tree Sequence)
# By passing a dictionary to `samples`, msprime generates haploid samples 
# matching the logic of your original msprime.Sample(pop, 0) list.
ts = msprime.sim_ancestry(
    samples={"pop0": diploid_sample_pop0, "pop1": diploid_sample_pop1},
    demography=demography,
    sequence_length=seq_length,
    recombination_rate=rho
)


# 3. Simulate Mutations
mutated_ts = msprime.sim_mutations(ts, rate=mu)


# 4. Write output to VCF
# Ploid=2 automatically pairs the generated haploid nodes into diploids 
# and assigns them the names generated in sampleIDs.
with sys.stdout as vcffile:
    mutated_ts.write_vcf(
        vcffile, 
        contig_id=chromID, 
        individual_names=sampleIDs
    )
