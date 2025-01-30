import numpy as np
from math import comb, e
from scipy.sparse import csr_matrix, identity
import scipy.sparse.linalg as sla
import scipy.linalg as la
import tskit as ts
from statsmodels.distributions.empirical_distribution import ECDF
from multiprocessing import Pool

def n_samps_mat(nlineages, n, ap, hp, r, mut):
# to avoid confusion, nlineages and n are both in terms of number of chromosomes. For a diploid organism, double the effective population size before inputting.
    stateslist = []
    dim1 = ((nlineages+1)*(nlineages+2)/2) - 3
    for i in range(nlineages+1):
        for j in range(nlineages+1-i):
            stateslist.append((i,j))
    stateslist = stateslist[1:]
#    print(stateslist)
    vals = []
    cols = []
    rows = []
    nnz = 0
    for state in stateslist:
        rows.append(nnz)
        for i in stateslist:
            first = True
            if state == (0,1) or state == (1,0):
                continue
            elif i == (state[0]-1, state[1]+1) and state[0] >= 1:
                vals.append(state[0]*r*ap/(state[0]+state[1]))
            elif i == (state[0]+1, state[1]-1) and state[1] >= 1: 
                vals.append(state[1]*r*ap/(state[0]+state[1]))
            elif i == (state[0] - 1, state[1]) and state[0] >= 2:
                vals.append((1/(hp*n))*comb(state[0], 2)/(state[0]+state[1]))
            elif i == (state[0], state[1]-1) and state[1] >= 2:
                vals.append((1/(hp*n))*comb(state[1], 2)/(state[0]+state[1]))
            else:
                continue
            nnz += 1
            cols.append(stateslist.index(i))
    rows.append(nnz)
    for i in range(len(rows)-1):
        if rows[i] == rows[i+1]:
            continue
        else:
            a =(-1 * sum(vals[rows[i]:rows[i+1]]))
            if i > max(cols[rows[i]:rows[i+1]]):
                vals.insert(rows[i+1], a)
                cols.insert(rows[i+1], i)
                for k in range(i+1,len(rows)):
                    rows[k] +=1
            else:
                for j in cols[rows[i]:rows[i+1]]:
                    if j > i:
                        ind = rows[i] + cols[rows[i]:rows[i+1]].index(j) 
                        vals.insert(ind, a)
                        cols.insert(ind, i)
                        for k in range(i+1,len(rows)):
                            rows[k] +=1
                        break
                    else:
                        continue
#    print(vals)
#    print(cols)
#    print(rows)
    rows = np.unique(rows)
    last = stateslist.index((1,0))
    #print(last)
    first = stateslist.index((0,1))
    #print(first)
#supposed to chop out coalescent states but doesn't work rn
    for i in range(len(rows)):
        if rows[i] > cols.index(first):
            rows[i] -= 1
        if rows[i] > cols.index(last):
            rows[i] -= 1
    del vals[cols.index(first)]
    for i in range(len(cols)):
        if cols[i] > first:
            cols[i] -= 1
    cols.remove(first)
    last -= 1
    ind = cols.index(last)
    del vals[ind]
    for i in range(len(cols)):
        if cols[i] > last:
            cols[i] -= 1
    del cols[ind]
    mat = csr_matrix((vals, cols, rows))
    #print(mat.toarray())
    #return mat
    return (identity(dim1, format="csr") - mat*(1/mut)).toarray()
def neut_mat(n, nlineages, mut):
    coal_rates = [-comb(x,2)/n for x in range(2,nlineages+1)]
    coal_rates.reverse()
    mat1 = np.diag(coal_rates) + -1*np.diag(coal_rates[:-1], k=1)
    l = mat1.shape[0]
    return (identity(l) - mat1*(1/mut))


def matmean(mat, vec):
    inv = la.inv(mat)
    ln = inv.shape[0]
    idd = np.eye(ln)
    vec2 = vec.dot(inv)
    vec3 = la.inv(idd - inv).dot(np.ones(ln))
    return vec2.dot(vec3)
def matmean_cont(mat, vec):
    inv = sla.inv(-1*mat).toarray()
    return vec.dot(inv.dot(np.ones(vec.shape[0])))

def cdf(mat, vec, x):
    mat = np.linalg.inv(mat)
    vec = vec.dot(mat)
    matpow = np.linalg.matrix_power(mat, x)
    vec = np.ravel(vec)
    ones = np.ones(vec.shape[0])
    return 1 - vec.dot(np.ravel(matpow.dot(ones)))
def distance_to_prob(distance, recomb_rate):
    distance = float(distance)
    x = 0.5*(1-e**(-2*distance*recomb_rate))
    return x
def prob_to_distance(prob, recomb_rate):
    prob = float(prob)
    # print(type(recomb_rate))
    x=(-(1/(2*recomb_rate))*np.log(1-2*prob))
    return x
def sampler(state, a0, a1):
    return np.concatenate((np.random.choice(a0, state[0], replace = False), np.random.choice(a1, state[1], replace = False)))

def pre_ecdf(treename, samp_size, num_samps, num_windows, recomb_rate, end, num_tasks):
    #need to redo windowing code to save window size, be less messy, and incorporate sensible endpoints
    #probably just do one round of windowing and then just add up windows that need to be paired together.
    tree = ts.load(str(treename))
    for x in tree.variants(left=int(1e6), right=int(1e6)+1):
        g = x
    alleles = np.unique(g.genotypes, return_counts=True)
    print(alleles)
    kept=[]
    tmp_count = 0
    kept = [alleles[0][x] for x in range(len(alleles[0])) if alleles[1][x] >= tree.samples().shape[0]*0.05]
    print(kept)
    counts = [alleles[1][np.where(alleles[0] == int(x))][0] for x in kept]
    tot = sum(counts)
    print(counts)
    freqs = [round(samp_size*x/tot) for x in counts]
    if len(kept) > 2:
        raise Exception("More than 2 non-trivial alleles")
    state = (freqs[0], freqs[1])
    print(state)
    a0 = [x for x in tree.samples() if g.alleles[g.genotypes[x]] == g.alleles[kept[0]]]
    print(len(a0))
    a1 = [x for x in tree.samples() if g.alleles[g.genotypes[x]] == g.alleles[kept[1]]]
    print(len(a1))
    prob_breakpoints = np.logspace(-5, np.log10(end), int(num_windows))
    avgs = [np.amin([prob_breakpoints[x], prob_breakpoints[x+1]]) for x in range(len(prob_breakpoints)-1)]
    #print(prob_breakpoints)
    distance_breakpoints = [prob_to_distance(x, recomb_rate) for x in prob_breakpoints]
    distance_breakpoints1 = [(x + 1e6) for x in distance_breakpoints]
    distance_breakpoints2 = [(1e6-x) for x in distance_breakpoints]
    distance_breakpoints1.insert(-1, tree.sequence_length)
    distance_breakpoints2.insert(0, 0)
    distance_breakpoints2.insert(-1, 1e6)
    distance_breakpoints3 = distance_breakpoints2 + distance_breakpoints1
    distance_breakpoints3.sort()
    window_sizes = [distance_breakpoints[x+1] - distance_breakpoints[x] for x in range(len(distance_breakpoints)-1)]
    print(distance_breakpoints3)
    #need to multi-thread this part
    with Pool(num_tasks) as pool:
        samps_list = pool.starmap(sampler, [(state, a0, a1) for i in range(num_samps)])
    #samps_list = np.array(samps_list)
    #samps_list = [np.concatenate((np.random.choice(a0, state[0], replace = False), np.random.choice(a1, state[1],replace=False))) for x in range(num_samps)]
    samp_sites = np.array([tree.segregating_sites(sample_sets = samp, span_normalise=False, windows=distance_breakpoints3) for samp in samps_list])
    winnum = int((len(distance_breakpoints3)-1)/2)
    samp_sites = [samp_sites[:,winnum-x] + samp_sites[:,winnum+x] for x in range(winnum-1)]
    #multi-thread till here
    stateslist = []
    for i in range(int(samp_size)+1):
        for j in range(int(samp_size)+1-i):
            stateslist.append((i,j))
    stateslist = stateslist[1:]
    stateslist.remove((1,0))
    stateslist.remove((0,1))
    index = stateslist.index(state)
    vec = np.zeros(len(stateslist))
    #print(state)
    #print(stateslist)
    vec[index] = 1
    #samp_sites = np.transpose(samp_sites)
    return [np.array(samp_sites), vec, window_sizes, avgs]
def ecdf_arr(samp_sites):
    tmp = []
    for i in range(samp_sites.shape[0]):
        ecdf = ECDF(samp_sites[i])
        list1 = [ecdf(x) for x in np.arange(np.amin(samp_sites[i]), np.amax(samp_sites[i])+1)]
        #store min val in ecdf for easy reconstruction of bounds
        list1.append(np.amin(samp_sites[i]))
        #include 0 for one less than min number of seg sites
        list1.insert(0,0)
        tmp.append(list1)
    return tmp
def ks_dist(mat, vec, ecdf_vec):
    mini = round(ecdf_vec[-1])
    ecdf_vec = ecdf_vec[:-1]
    #print(ecdf_vec)
    maxi = len(ecdf_vec) + mini - 1
    acdf_vec = np.array([cdf(mat, vec, x) for x in range(int(mini-1), int(maxi))])
    #print(acdf_vec)
    return np.amax(np.abs(np.subtract(acdf_vec, ecdf_vec)))
def ks_solve(ecdf_vec, win_size, vec, lins, n, ap, hp, r, mut):
    mat =  n_samps_mat(lins, n, ap, hp, r, mut*win_size)
    return ks_dist(mat, vec, ecdf_vec)
def neut_ks(ecdf_vec, win_size, lins, n, mut):
    mat = neut_mat(n, lins, mut*win_size)
    vec =  np.zeros(lins-1)
    vec[-1] = 1
    return ks_dist(mat, vec, ecdf_vec)
