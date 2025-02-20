import pyslim, tskit, msprime
import sys

file1 = str(sys.argv[1])
N = float(sys.argv[2])
mu = float(sys.argv[3])
r = float(sys.argv[4])
name = str(sys.argv[5])
ts = tskit.load(file1)
recap = pyslim.recapitate(ts, ancestral_Ne=N, recombination_rate=r, random_seed=5)
next_id = pyslim.next_slim_mutation_id(recap)
mts = msprime.sim_mutations(recap, rate=mu, keep=True, model=msprime.SLiMMutationModel(type=10, next_id=next_id))
mts.dump("./" + name + "_recap.trees")
