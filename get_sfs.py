from math import e
import tskit as ts
import numpy as np
def distance_to_prob(distance, recomb_rate):
    #helper function, implements Haldane mapping function
    distance = float(distance)
    x = 0.5*(1-e**(-2*distance*recomb_rate))
    return x
def prob_to_distance(prob, recomb_rate):
    #inverse of Haldane mapping function
    prob = float(prob)
    # print(type(recomb_rate))
    x=(-(1/(2*recomb_rate))*np.log(1-2*prob))
    return x
def max_dist(lst):
    #finds longest sequence run between two recombination breakpoints
    #I did get lazy and have AI write this, as you can see by it actually raising specific errors instead of generic exceptions
    if len(lst) < 2:
        raise ValueError("The list must contain at least two elements.")
    max_diff = 0
    elements = (None, None)
    for i in range(len(lst) - 1):
        current_diff = abs(lst[i] - lst[i + 1])
        if current_diff > max_diff:
            max_diff = current_diff
            elements = (lst[i], lst[i + 1])
    return elements

def find_loci(breakslist, pos):
    not_found = True
    low, high = 0, len(breakslist) - 1
    while not_found:
        mid = (low + high)//2
        if (low == (high - 1)):
            return (low, high)
        if breakslist[mid] == pos:
            return [mid]
        if breakslist[mid] < pos:
            low = mid
        else:
            high = mid

def get_mean_sfs_loc(treename, samp_size, num_samps, num_loci):
    tree = ts.load(str(treename))
    tree_points = tree.breakpoints(as_array=True)
    posind = find_loci(tree_points, 1e6)[0]
    locilist = [x - num_loci + posind for x in range(2*num_loci+2)]
    poslist = [0, *[tree_points[x] for x in locilist], 2e6+1]
    samps_list = np.random.rand(num_samps, len(tree.samples())).argpartition(samp_size,axis=1)[:,:samp_size]
    samp_sites = np.array([tree.allele_frequency_spectrum(sample_sets = [samps_list[i]], span_normalise=False, windows=poslist, polarised=True) for i in range(samps_list.shape[0])]) 
    return np.sum(samp_sites, axis=0)[1:-1, 1:-1]/num_samps

def get_mean_sfs_win(treename, samp_size,  num_samps, dist):
    tree = ts.load(str(treename))
    poslist = [0, 1e6 - dist, 1e6 + dist,  2e6+1]
    samps_list = np.random.rand(num_samps, len(tree.samples())).argpartition(samp_size,axis=1)[:,:samp_size]
    samp_sites = np.array([tree.allele_frequency_spectrum(sample_sets = [samps_list[i]], span_normalise=False, windows=poslist, polarised=True) for i in range(samps_list.shape[0])]) 
    return np.sum(samp_sites, axis=0)[1:-1, 1:-1]/num_samps

def get_mean_sfs_win_time_series(treename, n, samp_size,  num_samps, dist):
    tree = ts.load(str(treename))
    poslist = [0, 1e6 - dist, 1e6 + dist,  2e6+1]
    samps = tree.samples()
    samp_sets = np.reshape(samps, (int(len(samps)/n), n))
    indice_sets = [np.random.rand(num_samps, len(tree.samples())).argpartition(samp_size,axis=1)[:,:samp_size] for set in samp_sets]
    samps_list = [np.take(samp_sets[i], indice_sets[i]) for i in range(samp_sets.shape[0])]
    samp_sfss = [np.array([tree.allele_frequency_spectrum(sample_sets = [samp_list[i]], span_normalise=False, windows=poslist, polarised=True) for i in range(samp_list.shape[0])]) for samp_list in samps_list]
    return [np.sum(samp_sfs, axis=0)[1:-1, 1:-1]/num_samps for samp_sfs in samp_sfss]