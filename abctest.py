import scipy.stats as sit
import pyabc
import n_samps_model as n 
import numpy as np
import os
import tempfile
import sys
from scipy.optimize import brentq
treename = str(sys.argv[1])
nn = int(sys.argv[2])
lins = int(sys.argv[3])
pre_ecdf = n.pre_ecdf(treename, lins, 1000, 5, 1e-6, 1e-3)
arr = n.ecdf_arr(pre_ecdf[0])
obs = {str(i): 0 for i in range(len(pre_ecdf[3]))}

print(obs)
class Prior1(pyabc.DistributionBase):
    def __init__(self):
        self.ap = pyabc.RV("uniform", 0.06, 0.94)
        self.hp = pyabc.RV("uniform", 0, 0.95)
    def rvs(self, *args, **kwargs):
        while True:
                ap, hp = self.ap.rvs(), self.hp.rvs()
                if hp < (ap - 0.05):
                    return pyabc.Parameter(ap=ap, hp=hp)
    def pdf(self, x):
        ap, hp = x["ap"], x["hp"]
        if hp >= (ap -0.05):
                return 0.0
        return self.ap.pdf(ap) * self.hp.pdf(hp)
symrv = Prior1()
n_list = [n.neut_ks(arr[x], pre_ecdf[2][x], lins, nn, 1e-6) for x in range(len(pre_ecdf[3]))]
n_dict = {str(i): n_list[i] for i in range(len(n_list))}
def sym(dict1):
    ap = dict1["ap"]
    hp = dict1["hp"]
    list1 = [n.ks_solve(arr[x], pre_ecdf[2][x], pre_ecdf[1], lins, nn, ap, hp, pre_ecdf[3][x], 1e-6) for x in range(len(pre_ecdf[3]))]
    return {str(i): list1[i] for i in range(len(list1))}
def over(dict2):
    p = dict2["p"]
    fudge = dict2["fudge"]
    list1 = [n.ks_solve(arr[x], pre_ecdf[2][x], pre_ecdf[1], lins, nn, p, p-fudge, pre_ecdf[3][x], 1e-6) for x in range(len(pre_ecdf[3]))]
    return {str(i): list1[i] for i in range(len(list1))} 
def neut(dict3):
    return n_dict
def dist(dict1, dict2):
    a = np.array(list(dict1.values()))
    b = np.array(list(dict2.values()))
    return np.linalg.norm(a-b)
models = [sym, over, neut]
priors = [symrv, pyabc.Distribution(p=pyabc.RV("uniform", 0, 1), fudge=pyabc.RV("uniform", 0, 0.05)), pyabc.Distribution(k=pyabc.RV("uniform", 0, 1))]

abc = pyabc.ABCSMC(models, priors, pyabc.distance.PNormDistance(p=2))
db_path = "sqlite:///" + os.path.join(tempfile.gettempdir(), "test0.db")
bflist=[]
for i in range(3):
    while True:
        try:
            history = abc.new(db_path, obs)
            history = abc.run(minimum_epsilon=0.5, max_nr_populations=7)
            break
        except AssertionError:
            pass

    prob=(history.get_model_probabilities())
    print(prob)
    nonzero = False
    i = -1
    while nonzero==False:
        if prob.iloc[i][1] != 0:
            row = prob.iloc[i]
            bf = row[0]/row[1]
            nonzero = True
        else: 
            i -= 1
    print(bf)
    bflist.append(bf)
print(f'Log10 of Bayes Factor = {np.log10(np.mean(bflist))}')



