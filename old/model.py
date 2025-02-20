import numpy as np
from numpy import exp
from scipy.integrate import quad, dblquad
from scipy import linalg
from math import factorial
from scipy.optimize import fsolve
import sys, warnings
from multiprocessing import Pool
#from mpi4py.futures import MPIPoolExecutor as mpiPool

s1 =np.ones((3,1))
rng = np.random.default_rng()

def g_mat(n, p, r):
#generates the sub-intensity matrix for the Markov chain
    c = 1/(2*n*p)
    w = - c - r 
    return np.array([[w, r, 0], [0.5*r,  -1*r, 0.5*r], [0, r, w]])
def g_mat1(n, p, r):
    #attempt at using infinite series to improve model fit, did not work
    c = 1/(2*n*p)
    v1 = [0, r/(1-(r**2)),  (r**2)/(1-(r**2)), 1/(1-c)]
    v1[0] = -1*sum(v1)
    v1 = v1[:-1]
    v2 = [2*r/(4-r**2), 0, 2*r/(4-r**2)]
    v2[1] = -1*sum(v2)
    v3 = v1[::-1]
    return np.array([v1, v2, v3])

def state_vec(s,g):
    #computes HW-eq at harmonic mean allele freq
    p = (4+3*s)/(4+s)
    p0 = 1/(1+(p**(-1*g/2)))
    return np.array([[p0**2, 2*(1-p0)*p0, (1-p0)**2]])
def state_vec_p(p):
    #computes HW-eq at current allele freq
    return np.array([p**2, 2*(1-p)*p, (1-p)**2])

def ph_cdf(mat, state_vec1, x):
    #cdf of the coalescence time (denoted PH since this is a phase-type distribution)
    return 1 - state_vec1.dot(linalg.expm(x*mat).dot(s1))
def ph_pdf(mat, state_vec1,x):
    #pdf of coalescence time 
    if x < 0:
        return 0
    else:
        s = -1*mat.dot(s1)
        return state_vec1.dot(linalg.expm(x*mat).dot(s))
def ph_pdf_prime(mat, vec, x):
    #derivative of pdf of coalescence time
    s = -1*mat.dot(s1)
    return vec.dot(mat.dot((linalg.expm(x*mat).dot(s))))

def ex_t(mat, vec1):
    #expectation of coalescence time
    U =linalg.inv(-1*mat)
    u1= U.dot(s1)
    vec2 = vec1
    return (vec2.dot(u1))[0]

def ph_qf(mat, vec, t):
    tmp = lambda x: ph_cdf(mat, vec, x) - t
    return fsolve(tmp, ex_t(mat, vec))[0]


def ph_rv(mat, vec):
    #random sample of coalescence time
    return ph_qf(mat, vec, np.random.uniform())
def dist_rv(mat, vec, mu,n):
    #random sample of diversity
    t = ph_rv(mat, vec)
    return np.random.normal(2*mu*t, np.sqrt(4*mu*t/n))

def single_loc_pval(mat, vec, obs_div, mu, n, nruns):
    #two tailed p-val for the observed diversity assuming balancing selection at a locus with a recombination probability of r relative to our test locus
    #r is captured in the input matrix
    samps = np.sort([dist_rv(mat, vec, mu, n) for x in range(nruns)])
    std = np.std(samps)
    mean = np.mean(samps)
    distance = abs((obs_div - mean)/std)
    beginning, end = mean-(std*distance), mean+(std*distance)
    arr = samps[np.logical_or(samps < beginning, samps > end)]
    return arr.shape[0]/nruns

def bisect(f, a, b, tol):
    if f(a) * f(b) >= 0:
        print("error")
        return
    c = a
    while (b-a) >= tol:
        c = (a+b)/2
        if f(c) == 0.0:
            break
        elif f(c)*f(a) < 0:
            b = c
        else:
            a = c
    return c

def allele_freq(s, t, g):
    #models the allele frequency at a given t in the cycle, where g is the generations per season
    p = (4+3*s)/(4+s)
    i = t % (2*g)
    if i < g:
        return 1/(1+p**(i - 0.5*g))
    else:
        return 1/(1+p**(1.5*g - i))

def big_f(s, g, r):
    #this is a scaling factor for population size based on hitchhiking
    def delta_i(s_1, t_1, g_1):
        return (allele_freq(s_1, t_1 + 1, g_1) - allele_freq(s_1, t_1, g_1))
    def s_i(s, t, g, r):
        def summand_s(s, t, g, r, j):
            return (delta_i(s, t+j, g)*(exp(-r*j)/(1-exp(-r*2*g))))
        return sum(summand_s(s, t, g, r, x) for x in range(0, 2*g))
    def summand_t(s, t, g, r):
        pi_i = allele_freq(s, t, g)
        return ((1/(pi_i * (1 - pi_i)))*(s_i(s, t, g, r))**2)
    sum1 = sum(summand_t(s, x, g, r) for x in range(0, 2*g))
    f = 1 + sum1/(2*g)
    return f
def new_big_f(freqarray, delta_array, g, r):
    #g = int(g)
    #print(type(g))
    #print(g)
    def s_i(delta_array, t, g, r):
        def summand_s(delta_array, t, g, r, j):
            ind = (t+j)%(2*g)
            return (delta_array[ind]*(exp(-r*j)/(1-exp(-r*2*g))))
        return sum((summand_s(delta_array, t, g, r, x) for x in range(2*g)))
    #print([s_i(delta_array, t, g, r) for t in range(2*g)])
    def summand_t(freqarray, delta_array, t, g, r):
        pi_i = freqarray[t]
        return ((1/(pi_i * (1 - pi_i)))*(s_i(delta_array, t, g, r))**2)
   # print(type(g))
   # print(f'g: {g}')
    sums = [summand_t(freqarray, delta_array, x, g, r) for x in range(int(2*g))]
    #print(sums)
    sum1 = sum(sums)
    f = 1 + sum1/(2*g)
    return f

def avg_cdf(r0, r1, n, p, x, vec):
    q = lambda j: ph_cdf(g_mat(n, p, j), vec, x)
    return quad(q, a=r0, b=r1)[0]/(r1-r0)
def avg_cdf_bio_explicit(r0, r1, n, x, vec, s, g):
    q = lambda j: ph_cdf(g_mat(n, 1/(2*big_f(s,g,j)), j), vec, x)
    return quad(q, a=r0, b=r1)[0]/(r1-r0)
def band_qf(vec, r0, r1, n, s, g, t):
    tmp = lambda x: avg_cdf_bio_explicit(r0, r1, n, x, vec, s, g) - t
    return fsolve(tmp, 2*n)[0]


def band_rv(vec, r0, r1, n, s, g):
#   random sample of coalescence time
        return band_qf(vec, r0, r1, n, s, g, np.random.uniform())
     
def band_dist_rv(vec, r0, r1,s, g, mu,n):
    #random sample of diversity
    t = band_rv(vec, r0, r1, n, s, g)
    return np.random.normal(2*mu*t, np.sqrt(4*mu*t/n))


def avg_cdf_bio_explicit2(r0, r1, n, x, vec, freqarray,delta_array, g):
    q = lambda j: ph_cdf(g_mat(n, 1/(2*new_big_f(freqarray, delta_array,g,j)), j), vec, x)
    return quad(q, a=r0, b=r1)[0]/(r1-r0)

def band_qf2(vec, r0, r1, n, freqarray, delta_array, g, t):
    tmp = lambda x: avg_cdf_bio_explicit2(r0, r1, n, x, vec, freqarray, delta_array, g) - t
    return fsolve(tmp, 2*n)[0]


def band_rv2(vec, r0, r1, n, freqarray, delta_array, g):
    #random sample of coalescence time
    return band_qf2(vec, r0, r1, n, freqarray, delta_array, g, np.random.uniform())

def band_dist_rv3(vec, r0, r1,freqarray, delta_array, g, mu,n):
    #random sample of diversity
    t = band_rv2(vec, r0, r1, n, freqarray, delta_array, g)
    return np.random.normal(2*mu*t, np.sqrt(4*mu*t/n))

def band_dist_rv_2(vec, r0, r1, s, g, mu, n, runs=1000):
    end, beginning = np.log10(r0), np.log10(r1)
    points = np.logspace(end, beginning, int(runs))
    mats = [g_mat(n, 1/(2*big_f(s,g,j)), j) for j in points]
    samps = [dist_rv(mat, vec, mu, n) for mat in mats]
    return np.mean(samps)

def subintmat_preparer(s,g,r,n):
    c = big_f(s,g,r)/n
    val = np.sqrt(c**2 + 4*(r**2))
    ev1, ev2, ev3 = (-c -r), (-val -c - 2*r)/2, (val -c - 2*r)/2
    #exps = [np.exp(y) for y in [ev1, ev2, ev3]]
    mat1 = np.array([[-1,1,1], [0, (c - val)/(2*r), (c + val)/(2*r)], [1,1,1]])
    mat2 = np.array([[-0.5, 0, 0.5], [c/(4*val) + 0.25, -1*r/val, c/(4*val) + 0.25], [-c/(4*val) + 0.25, r/val, -c/(4*val) + 0.25]])    
    #diag = np.diag(exps)
    return (mat1, [ev1, ev2, ev3], mat2)
def subintmat_exp(matslist, x):
     diag=np.diag([exp(y*x) for y in matslist[1]])
     return matslist[0].dot(diag.dot(matslist[2]))
def subintmat_exp2(s,g,r,n, x):
    c = big_f(s,g,r)/n
    val = np.sqrt(c**2 + 4*(r**2))
    y1, y2, y3 = exp(x*(-c -r)), exp(x*(-val -c - 2*r)/2), exp(x*(val -c - 2*r)/2)
    a = c*y2 - c*y3 + val*y2 + val*y3
    diag = 2*a + 4*val*y1
    antidiag = diag - 8*val*y1
    vert = -8*(r*y2 - r*y3)
    sides = ((c**2)*y2 - (c**2)*y3 - (val**2)*y2 + (val**2)*y3)/r
    return (1/(8*val))*np.array([[diag, vert, antidiag], [sides, -4*(c*y2-c*y3-val*y2-val*y3), sides], [antidiag, vert, diag]]) 



def new_ph_cdf(matslist, vec, x):
    return 1 - vec.dot(subintmat_exp(matslist,x).dot(s1))
#def new_dist_rv

from time import time
def mat_exp_test(runs):
    tlist1 = []
    tlist2 = []
    #difflist = []
    matslist = subintmat_preparer(0.5, 20, 1e-5, 1e4)
    vec = state_vec_p(0.5)
    for i in range(runs):
        t1 = time()
        p = 1/(2*big_f(0.5, 20, 1e-5))
        mat = g_mat(1e4, p, 1e-5)
        mat = mat*i
        linalg.expm(mat)
        t2 = time() - t1
        tlist1.append(t2)
        t3 = time()
        subintmat_exp2(0.5, 20, 1e-5, 1e4, i)
        t4 = time() - t3
        tlist2.append(t4)
        #difflist.append(a-b)
    return (sum(tlist1)/sum(tlist2))


#these approximations are all for the symmetric case right now, but I will extend them to the general case soon

def hmaf(s,g):
    #harmonic mean allele frequency
	p=(4+3*s)/(4+s)
	return (2*g*(p-1))/(2*g*(p-1)+(p**(g/2)-p**(-g/2))*(p+1))

def change_over_n_het(r, s,g,p, n, u):
    #does change over neutral heterozygosity, I don't really use this as diversity through Nei's pi is easier to work with
	f = big_f(s, g, r)
	h20 = 1-(1/((8*n*u/f)+1))
	h_11 = 1-(1/(4*u*((1/r) + (2*n/f))+1))
	het = (p**2 + (1-p)**2)*h20 + p*(1-p)*h_11 
	nhet = (4*n*u)/(4*n*u +0.5)
	return het/nhet



def change_over_n_div(r, s, g, n, p1):
   #returns ratio of (expected) diversity over neutral expectation, r is the recombination probability with the focal locus
    #if (r <= 1/n):
     #   p = hmaf(s,g)
    #else:
    p = 1/(2*big_f(s,g,r))
    mat = g_mat(n, p, r)
    vec1 = state_vec_p(p1)
    t = ex_t(mat, vec1)
    return t/(2*n)
def one_loc_pval(r,s,g,n,u,p,nruns, obs_div):
	#takes input of diversity, not ratio
	#just a function wrapper, but runs directly from biological parameters
	pm = 1/(2*big_f(s,g,r))
	mat = g_mat(n, pm, r)
	vec = state_vec_p(p)
	return single_loc_pval(mat, vec, obs_div, u, n, nruns)



def avg_change_div_num(s1, g1, r1, r2, n1, p1):
    #just window average diversity by integration
	a = quad(change_over_n_div, a=r1, b=r2, args=(s1,g1,n1, p1) )
	dist = r2 - r1
	return [a[0]/dist, a[1]/(dist**0.5)]
def partitioner(whole, n):
    whole, n = int(whole), int(n)
    rem = whole % n
    a = int((whole - rem)/n)
    list1 = [a for x in range(n)]
    for i in range(rem):
        list1[i] += 1
    return list1
    

def many_rvs(vec, windows,freqarray, delta_array, g, u,n, x):
    np.random.seed()
   
    '''
    print(f"s:{s}")
    print(f"g:{g}")
    print(f"u:{u}")
    print(f"n:{n}")
    print(f"x:{x}")
    '''
    return [[band_dist_rv3(vec, j[0], j[1], freqarray, delta_array, g, u, n) for i in range(x)] for j in windows]

def whole_chrom_test(w_vec, d_vec, s, g, n, u, p, nruns, ntasks, strapfac=1, mode="div"):
    #w_vec - list of window breakpoints (recombination probability with focal locus)
    #d_vec - list of average diversities in each window
    rng = np.random.default_rng()
    if mode=="ratio":
        d_vec = [x*4*n*u for x in d_vec]
    freqarray = [allele_freq(s, x, g) for x in range(2*g)]
    delta_array = [freqarray[(i+1)%(2*g)] - freqarray[i%(2*g)] for i in range(2*g)]
    windows = [(w_vec[i], w_vec[i+1]) for i in range(len(w_vec)-1)]
    #print(d_vec)
    #print(windows)
    ex_div = [avg_change_div_num(s, g, x[0], x[1], n, p)[0]*4*n*u for x in windows]
    #print(ex_div)
    vec = state_vec_p(p)
    partitions = partitioner(nruns, ntasks) 
    with Pool(ntasks) as pool:
        results = pool.starmap(many_rvs, [(vec, windows, freqarray, delta_array, g, u, n, partitions[i]) for i in range(len(partitions))])
    div_samps = np.array(results)
    div_samps = [np.ravel(div_samps[:,i]) for i in range(len(windows))]
    #div_samps = [[band_dist_rv(vec, x[0], x[1], s, g, u, n) for i in range(nruns)] for x in windows
    #print(np.mean(div_samps))
    div_vecs = [[rng.choice(x) for x in div_samps] for i in range(strapfac*nruns)]
    #print(div_vecs)
    m_dist = np.sort([np.mean([abs(y) for y in np.subtract(x, ex_div)]) for x in div_vecs])
    #print(m_dist)
    np.savetxt("m_samps.txt",m_dist)
    observed_m = np.mean([abs(x) for x in np.subtract(ex_div, d_vec)])
    print(f"observed error: {observed_m}")
    std = np.std(m_dist)
    mean = np.mean(m_dist)
    distance = np.abs((observed_m - mean)/std)
    print(f'distance: {distance}')
    #beginning, end = mean-(std*distance), mean+(std*distance)
    arr = m_dist[m_dist > observed_m]
    return arr.shape[0]/(nruns*strapfac)

if __name__ == "__main__":
    warnings.simplefilter("ignore")
    points = np.loadtxt("./WFrun1/breakpoints_" + str(sys.argv[5]) + ".gz")
    div_vec = [x*4e-2 for x in np.loadtxt("./WFrun1/diversities_" + str(sys.argv[5]) + ".gz")]
    print("loaded")
    a=whole_chrom_test(points, div_vec, 0.5, int(sys.argv[4]), 1e4, 1e-6, float(sys.argv[3]), int(sys.argv[1]), int(sys.argv[2]), strapfac=int(sys.argv[6]))
    print(f'p-val = {a}')
