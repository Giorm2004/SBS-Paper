import numpy as np 
from math import factorial as fac
from scipy.stats import poisson, binom, poisson_binom, Normal
from math import e 
from fast_poibin import PoiBin

def one_locality_error(D, r, cutoff):
    #maybe vectorize?
    return 1 - sum([((1-r)**x  +x*r*(1-r)**(x-1))*poisson.pmf(x, D) for x in range(2, cutoff)])

print(one_locality_error(200, 0.001, 800))
print(1 - binom.cdf(1, 200, 0.001))

def many_locality_error(Dlist, r, n_obs):
   l = [one_locality_error(D, r, 2*D) for D in Dlist]
   print(l)
   return poisson_binom.cdf(n_obs, l)


   
def final_binom_approx(k,n,p):
    if (n > 9*((1-p)/p)) and (n > 9*(p/(1-p))):
        return 1-Normal.cdf(k, n*p, np.sqrt(n*p*(1-p)))
    else:
        return 1 - poisson.cdf(k, n*p)





