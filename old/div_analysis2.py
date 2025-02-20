import pyslim, msprime, tskit, sys
import numpy as np
import matplotlib.pyplot as plt
from math import e
#gts = tskit.load("WFsim_" +sys.argv[1] + ".trees")
recomb_rate = float(sys.argv[2]) # per generation
#print(type(recomb_rate))
Ne = 10000 # generations
mut_rate = 1e-6 # per generation
#rts = pyslim.recapitate(gts,recombination_rate=recomb_rate, ancestral_Ne=Ne,random_seed=5)
#next_id = pyslim.next_slim_mutation_id(rts)
#rts = msprime.sim_mutations(rts, rate=mut_rate, model=msprime.SLiMMutationModel(type=10, next_id=next_id), keep=True,)
    
#rts.dump("WFsim_" + sys.argv[1] + ".trees")
rts = tskit.load("WFsim_" + sys.argv[1] + ".trees")
#prob_breakpoints = [1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.25]
def distance_to_prob(distance):
    distance = float(distance)
    x = 0.5*(1-e**(-2*distance*recomb_rate))
    return x
end = np.log10(distance_to_prob(rts.sequence_length-1e6))
breakpoints = np.linspace(-5, end, 15)
prob_breakpoints = [10**x for x in breakpoints]
print(prob_breakpoints)
def prob_to_distance(prob):
    prob = float(prob)
   # print(type(recomb_rate))
    x=(-(1/(2*recomb_rate))*np.log(1-2*prob))
    return x
distance_breakpoints = [prob_to_distance(x) for x in prob_breakpoints]
distance_breakpoints1 = [(x + 1e6) for x in distance_breakpoints]
distance_breakpoints2 = [(1e6-x) for x in distance_breakpoints]
distance_breakpoints1.insert(0, 0)
distance_breakpoints1.insert(1, 1e6)
distance_breakpoints1[-1] = rts.sequence_length
distance_breakpoints2.insert(0, 0)
distance_breakpoints2.insert(1, 1e6)
distance_breakpoints2[-1] = rts.sequence_length
distance_breakpoints2.sort()
divarr1 = rts.diversity(windows=distance_breakpoints1)
scaled_divarr1 = list(map(lambda x: x/(8*Ne*mut_rate), divarr1))[1:]
divarr2 = rts.diversity(windows=distance_breakpoints2)
scaled_divarr2 = list(map(lambda x: x/(8*Ne*mut_rate), divarr2))[:0:-1]
scaled_divarr = np.add(scaled_divarr1, scaled_divarr2).tolist()
np.savetxt("diversities_" + sys.argv[1] + ".gz", scaled_divarr)
np.savetxt("breakpoints_" + sys.argv[1] + ".gz", prob_breakpoints)
'''
graph_points = []
for i in range(len(breakpoints)):
    if (i+1 < len(breakpoints)):
        graph_points.append(0.5*(breakpoints[i]+breakpoints[i+1]))
plt.plot(graph_points, scaled_divarr)
plt.savefig("diversities_"+ sys.argv[1] + ".png")
'''


