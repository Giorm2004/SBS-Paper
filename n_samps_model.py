import numpy as np
from math import comb, e
from scipy.sparse import csr_matrix, identity
import scipy.sparse.linalg as sla
import scipy.linalg as la
import tskit as ts
from statsmodels.distributions.empirical_distribution import ECDF

def n_samps_mat(nlineages, n, ap, hp, r, mut, c=None):
# this generates the subintensity matrix for the PH distribution of coalescence times
    if c==None:
        c = 2*n*hp
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
                vals.append(state[0]*r*ap)
            elif i == (state[0]+1, state[1]-1) and state[1] >= 1: 
                vals.append(state[1]*r*ap)
            elif i == (state[0] - 1, state[1]) and state[0] >= 2:
                vals.append((1/(c))*comb(state[0], 2))
            elif i == (state[0], state[1]-1) and state[1] >= 2:
                vals.append((1/(c))*comb(state[1], 2))
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
    return mat
    #print(mat.toarray())
    #return mat
    #mat1 = (identity(dim1, format="csr") - mat*(1/(n*mut))).toarray()
    #return np.linalg.inv(mat1)
def s_transform(dim, lins, mat, mut):
    mat = mat.toarray()
    stateslist = []
    for i in range(int(lins)+1):
        for j in range(int(lins)+1-i):
            stateslist.append((i,j))
    stateslist = stateslist[1:]
    stateslist.remove((1,0))
    stateslist.remove((0,1))
    diag = np.diag(np.array([1/(state[0]+state[1]) for state in stateslist]))
    mat1 = (identity(dim) - mat*(1/(mut)))
    return np.linalg.inv(mat1)
def overdom_mat(nlineages, n, p, r):
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
                vals.append(state[0]*r*(1-p)/(2))
            elif i == (state[0]+1, state[1]-1) and state[1] >= 1:
                vals.append(state[1]*r*p/(2))
            elif i == (state[0] - 1, state[1]) and state[0] >= 2:
                vals.append((1/(2*p*n))*comb(state[0], 2))
            elif i == (state[0], state[1]-1) and state[1] >= 2:
                vals.append((1/(2*(1-p)*n))*comb(state[1], 2))
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
    return(mat)


def neut_mat(n, nlineages, mut):
    #matrix for regular  old coalescent, just reformulated as a PH distribution
    coal_rates = [-comb(x,2)/x for x in range(2,nlineages+1)]
    coal_rates.reverse()
    mat1 = np.diag(coal_rates) + -1*np.diag(coal_rates[:-1], k=1)
    l = mat1.shape[0]
    return np.linalg.inv((identity(l) - mat1*(1/(n*mut))))


def matmean(mat, vec):
    #computes mean of a DPH dist
    #inv = la.inv(mat)
    ln = mat.shape[0]
    idd = np.eye(ln)
    vec2 = np.linalg.inv(idd - mat).dot(np.ones(ln))
    return vec.dot(vec2)
def matmean_cont(mat, vec):
    #Mean of PH
    inv = sla.inv(-1*mat).toarray()
    return vec.dot(inv.dot(np.ones(vec.shape[0])))

def cdf(mat, vec, x):
    #cdf of DPH, need to rewrite in S+1 form as in Holboth
    #mat = np.linalg.inv(mat)
    #print(mat)
    matpow = np.linalg.matrix_power(mat, x+1)
    vec = np.ravel(vec)
    ones = np.ones(vec.shape[0])
    return 1 - vec.dot(np.ravel(matpow.dot(ones)))
def pmf(mat, vec, x):
    #pmf of DPH, need to rewrite in S+1 as in Holboth
    x = int(x)
    if isinstance(vec, list):
        vec = np.array(vec)
    s1 = np.ones(vec.shape[0])
    '''
    if mat.shape == (1,1):
        mat = mat[0,0]
        s2 = s1 - mat*s1
        return vec.dot(((mat**x)*s2))
    '''
    s2 = np.ravel(s1 - mat.dot(s1))
    matpow = np.linalg.matrix_power(mat, x)
    return vec.dot(np.ravel(matpow.dot(s2)))


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
def pre_ecdf(treename, samp_size, num_samps, num_windows, recomb_rate, end):
    #This processes a bunch of data that we need 
    tree = ts.load(str(treename))
    #Sometimes, we get an extra allele of a given type in SLiM, which throws off Allele indexing
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
    freq = counts[0]/tot
    if len(kept) > 2:
        raise Exception("More than 2 non-trivial alleles")
    #need to basically get rid of this. This was all about sampling so that the samples have an allele frequency as close as possible to that of the population, but this is just dumb
    #Instead, just choose samples randomly, and then just get state vector by using binomial distribution. 
    #By law of large numbers, this should just converge to whatever the true state vector is anyways. Basically generalized HW btw. 
    #first round of windowing, artificial
    #have to rework this so windows are far enough apart
    prob_breakpoints = np.logspace(-5, np.log10(end), int(num_windows))
    #print(prob_breakpoints)
    tree_points = tree.breakpoints(as_array=True)
    distance_breakpoints = [prob_to_distance(x, recomb_rate) for x in prob_breakpoints]
    print(distance_breakpoints)
    distance_breakpoints1 = [(x + 1e6) for x in distance_breakpoints]
    distance_breakpoints2 = [(1e6-x) for x in distance_breakpoints]
    distance_breakpoints1.insert(-1, tree.sequence_length)
    distance_breakpoints2.insert(0, 0)
    distance_breakpoints2.insert(-1, 1e6)
    distance_breakpoints3 = distance_breakpoints2 + distance_breakpoints1
    distance_breakpoints3.sort()
    window_sizes = [distance_breakpoints[x+1] - distance_breakpoints[x] for x in range(len(distance_breakpoints)-1)]
    print(distance_breakpoints3)
    #now we check all windows that recombination creates
    #if we try to estimate the coalescent of some locii with a recomb breakpoint inside it, we would need to take a convolution or something, which would explode our state space
    possible_windows = [[p for p in tree_points if distance_breakpoints3[i] <= p <= distance_breakpoints3[i+1]] for i in range(1,len(distance_breakpoints3)-2)]
    #if only one breakpoint in a window, we just arbitrarily chop off the sequence run at wherever the window ends:
    for l in possible_windows:
        #this is not actually good. Sometimes this breaks if the recomb breakpoint is right next to the artificial window end. What this should really do is:
        #if one breakpoint in window, add in whatever end of window maximizes length of seq run
        if len(l) == 1:
            a = 0
            for i in distance_breakpoints3:
                if l[0] < i-2:
                    a = int(i-1)
                    break
            l.append(a)
    print(possible_windows)
    max_dists = [max_dist(x) for x in possible_windows]
    #store the length of our new windows, need this to upscale mutation rate
    win_sizes = [tup[1] - tup[0] for tup in max_dists]
    #just gives us where on the genome these occur, need recomb probs to run anything
    recombdists = [distance_to_prob(abs(tup[1] - 1e6),recomb_rate) if tup[1] < 1e6 else distance_to_prob(tup[0] - 1e6, recomb_rate) for tup in max_dists]
    new_breaks = [item for tup in max_dists for item in tup]
    new_breaks.insert(0,0)
    new_breaks.insert(-1, tree.sequence_length)
    new_breaks.sort() 
    print(new_breaks)
    print("done with finding windows")
    samps_list = np.random.rand(num_samps, len(tree.samples())).argpartition(samp_size,axis=1)[:,:samp_size]
    print('done with samps_list')
    #samps_list = [np.concatenate((np.random.choice(a0, state[0], replace = False), np.random.choice(a1, state[1],replace=False))) for x in range(num_samps)]
    samp_sites = np.array(tree.segregating_sites(sample_sets = samps_list, span_normalise=False, windows=new_breaks)).T
    #winnum = int((len(distance_breakpoints3)-1)/2)
    samp_sites = np.array([samp_sites[:, x] for x in range(1,len(new_breaks)-1, 2)])
    samp_sites = [x.tolist() for x in samp_sites]
    #this just ensures we generate our state vector so states in matrix and vector match
    stateslist = []
    for i in range(int(samp_size)+1):
        for j in range(int(samp_size)+1-i):
            stateslist.append((i,j))
    stateslist = stateslist[1:]
    stateslist.remove((1,0))
    stateslist.remove((0,1))
    vec = np.zeros(len(stateslist)).tolist()
    #print(state)
    #print(stateslist)
    #do not do this. use binom instead
    for state in stateslist:
        if state[0]+state[1] == samp_size:
            vec[stateslist.index(state)] = comb(samp_size, state[0])*(freq**state[0])*((1-freq)**state[1])
    nvec = np.zeros(samp_size-1).tolist()
    nvec[-1] = 1
    #samp_sites = np.transpose(samp_sites)
    pre = [samp_sites, vec, win_sizes, recombdists, freq, nvec]
#    for i in pre:
 #       print(type(i))
  #      for j in i:
   #         print(type(j))
    return pre
def ecdf_arr(samp_sites):
    #returns ecdfs of segregating sites in each window
    tmp = []
    samp_sites = [np.array(x) for x in samp_sites]
    samp_sites = np.array(samp_sites)
    for i in range(samp_sites.shape[0]):
        ecdf = ECDF(samp_sites[i])
        list1 = [ecdf(x) for x in np.arange(np.amin(samp_sites[i]), np.amax(samp_sites[i])+1)]
        #store min val in ecdf for easy reconstruction of bounds
        list1.append(np.amin(samp_sites[i]))
        #include 0 for one less than min number of seg sites
        list1.insert(0,0)
        tmp.append(list1)
    return tmp
def eprob_arr(samp_sites):
    #returns epmfs of segregating sites in each window
    #need to actually allow 0 segregating sites here, which is solved by using Holboth S+1 form
    samp_sites = np.array(samp_sites)
    records = []
    for i in samp_sites:
        l = 1/(i.shape[0])
        unique, counts = np.unique(i, return_counts=True)
        counts = counts.astype(float)
        counts *= l
        records.append(np.array([unique, counts]))
    return records 

def kl_div(mat, vec, eprob_vec, siteslist):
    #1d KL div
    #should prob make multidimensional
    apmf_vec = np.array([pmf(mat, vec, x) for x in siteslist])
    return np.sum(np.array([eprob_vec[x] * np.log(eprob_vec[x]/apmf_vec[x]) for x in range(len(siteslist))]))
def ks_dist(mat, vec, ecdf_vec):
    #KS dist. Only well defined for 1d anyways. 
    mini = round(ecdf_vec[-1])
    ecdf_vec = np.array(ecdf_vec[:-1])
    print(ecdf_vec)
    maxi = len(ecdf_vec) + mini - 1
    acdf_vec = np.array([cdf(mat, vec, x) for x in range(int(mini-1), int(maxi))])
    print(acdf_vec)
    return np.amax(np.abs(np.subtract(acdf_vec, ecdf_vec)))
def ks_solve(ecdf_vec, win_size, vec, lins, n, ap, hp, r, mut):
    #this just computes the ks distance from parameter values directly, generates the mat
    #print(win_size)
    mat =  n_samps_mat(lins, n, ap, hp, r, mut*win_size)
    mat = s_transform(vec.shape, mat, mut*win_size)
    return ks_dist(mat, vec, ecdf_vec)
def sbs_kl_solve(eprob_tup, win_size, vec, lins, n, ap, hp, r, mut, c=None):
    #same as above but for KL div
    mat =  n_samps_mat(lins, n, ap, hp, r, mut*win_size, c)
    mat = s_transform(len(vec), lins, mat, mut*win_size)
    return kl_div(mat, vec, eprob_tup[1], eprob_tup[0])
def overdom_kl_solve(eprob_tup, win_size, vec, lins, n, p, r, mut):
    mat = overdom_mat(lins, n, p, r)
    mat = s_transform(len(vec), lins, mat, mut*win_size)
    return kl_div(mat, vec, eprob_tup[1], eprob_tup[0])
def neut_ks(ecdf_vec, win_size, lins, n, mut):
    #neutral
    mat = neut_mat(n, lins, mut*win_size)
    vec =  np.zeros(lins-1)
    vec[-1] = 1
    return ks_dist(mat, vec, ecdf_vec)
def neut_kl(tup, nvec, win_size, lins, n, mut):
    mat = neut_mat(n, lins, mut*win_size)
    nvec = np.array(nvec)
    return kl_div(mat, nvec, tup[1], tup[0])


#need to write neut_kl and also use binom state vector instead of unit state vector.
