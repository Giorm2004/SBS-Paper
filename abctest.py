import scipy.stats as sit
import pyabc
import n_samps_model as n 
import numpy as np
import os
import tempfile
import sys
#from scipy.optimize import brentq
treename = str(sys.argv[1])
nn = int(sys.argv[2])
lins = int(sys.argv[3])
db_num = str(sys.argv[4])
num_tasks = int(sys.argv[5])
pre_ecdf = n.pre_ecdf(treename, lins, 5000, 5, 1e-6, 1e-3, num_tasks)
arr = n.ecdf_arr(pre_ecdf[0])
obs = {str(i): 0 for i in range(len(pre_ecdf[3]))}

print(obs)
class Prior1(pyabc.DistributionBase):
    def __init__(self):
        self.ap = pyabc.RV("uniform", 0.06, 0.94)
        self.hp = pyabc.RV("uniform", 0.05, 0.95)
        self.ne = pyabc.RV("uniform", 0.75, 0.5)
    def rvs(self, *args, **kwargs):
        while True:
                ap, hp, ne = self.ap.rvs(), self.hp.rvs(), self.ne.rvs()
                if hp < (ap - 0.05):
                    return pyabc.Parameter(ap=ap, hp=hp, ne=ne)
    def pdf(self, x):
        ap, hp, ne = x["ap"], x["hp"], x["ne"]
        if hp >= (ap -0.05):
                return 0.0
        return self.ap.pdf(ap) * self.hp.pdf(hp) * self.ne.pdf(ne)
symrv = Prior1()
#n_list = [n.neut_ks(arr[x], pre_ecdf[2][x], lins, nn, 1e-6) for x in range(len(pre_ecdf[3]))]
#n_dict = {str(i): n_list[i] for i in range(len(n_list))}
def sym(dict1):
    ap = dict1["ap"]
    hp = dict1["hp"]
    ne = dict1["ne"]
    list1 = [n.ks_solve(arr[x], pre_ecdf[2][x], pre_ecdf[1], lins, ne*nn, ap, hp, pre_ecdf[3][x], 1e-6) for x in range(len(pre_ecdf[3]))]
    return {str(i): list1[i] for i in range(len(list1))}
def over(dict2):
    p = dict2["p"]
    fudge = dict2["fudge"]
    ne = dict2["ne"]
    list1 = [n.ks_solve(arr[x], pre_ecdf[2][x], pre_ecdf[1], lins, ne*nn, p, p-fudge, pre_ecdf[3][x], 1e-6) for x in range(len(pre_ecdf[3]))]
    return {str(i): list1[i] for i in range(len(list1))} 
def neut(dict3):
    ne = dict3["ne"]
    n_list = [n.neut_ks(arr[x], pre_ecdf[2][x], lins, ne*nn, 1e-6) for x in range(len(pre_ecdf[3]))]
    return {str(i): n_list[i] for i in range(len(n_list))}
def dist(dict1, dict2):
    a = np.array(list(dict1.values()))
    b = np.array(list(dict2.values()))
    return np.linalg.norm(a-b)
models = [sym, over, neut]
priors = [symrv, pyabc.Distribution(p=pyabc.RV("uniform", 0, 1), fudge=pyabc.RV("uniform", 0, 0.04), ne=pyabc.RV("uniform", 0.75, 0.5)), pyabc.Distribution(ne=pyabc.RV("uniform", 0.75, 0.5))]
print(len(priors))
abc = pyabc.ABCSMC(models, priors, pyabc.distance.PNormDistance(p=2))
db_path = "sqlite:///" + os.path.join(tempfile.gettempdir(), "test" + db_num + ".db")
bflist=[]
for i in range(3):
    while True:
        try:
            history = abc.new(db_path, obs)
            history = abc.run(minimum_epsilon=0.4, max_nr_populations=12)
            break
        except AssertionError:
            pass

    prob=(history.get_model_probabilities())
    print(prob)
    try:
        print(f'Neutral Prob: {prob.iloc[-1][2]}')
    except IndexError:
        raise Exception("Neutral model did not run")
    
    '''
    nonzero = False
    i = -1
    model_names = ["SBS", "Overdom", "Neutral"]
    for i in range(3):
        if prob.iloc[-1][i] == 1 or prob.iloc[-1][i] == 0:
            if prob.iloc[-1][i]==1:
                print(f'Model {model_names[i]} has reached a probability of 1.0!')
                for j in range(3):
                    if i != j:
                        nonzero = False
                        k = [-1]
                        while nonzero == False:
                            if prob.iloc[k][j] != 0:
                                bf = prob.iloc[k][i]/prob.iloc[k][j]
                                print(f'Last finite Bayes Factor for Prob {model_names[i]} over {model_names[j]} ={bf}'
                                nonzero == True
                            else:
                                k -= 1
                    bflist.append((bf, i, j))
            if prob.iloc[-1][i] == 0: 
                print(f'Model {model_names[i]} has reached a probability of 0.0!')
            else: 
                bf1, bf2 = (prob.iloc[-1][i]/prob.iloc[-1][(i-1)%3], i, (i-1)%3) , (prob.iloc[-1][i]/prob.iloc[-1][(i+1)%3], i, (i+1)%3)
                bflist.append(bf1)
                bflist.append(bf2)
    print(bf)
    bflist.append(bf)
print(f'Log10 of Bayes Factor = {np.log10(np.mean(bflist))}')
'''


