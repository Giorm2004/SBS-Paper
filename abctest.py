import scipy.stats as st
import pyabc
import n_samps_model as n 
import numpy as np
import os, json
import tempfile
import sys
#from scipy.optimize import brentq
treename = str(sys.argv[1])
nn = int(sys.argv[2])
lins = int(sys.argv[3])
db_num = str(sys.argv[4])
mut = float(sys.argv[5])
#pre_ecdf = n.pre_ecdf(treename, lins, 5000, 5, 1e-6, 1e-3, num_tasks)
#arr = n.ecdf_arr(pre_ecdf[0])
if os.path.exists(db_num+"_pre_ecdf.json"):
    with open(db_num + "_pre_ecdf.json", 'r') as f:
        pre_ecdf  = json.load(f)
else:
    pre_ecdf = n.pre_ecdf(treename, lins, 100000, 4, 1e-7, 1e-4)
    with open(db_num + "_pre_ecdf.json", 'w') as f:
        json.dump(pre_ecdf, f)
arr = n.eprob_arr(pre_ecdf[0])
vec1 = np.array(pre_ecdf[1])
#obs = {str(i): 0 for i in range(len(pre_ecdf[3]))}
obs = {'0': 0, '1': 0}
p = pre_ecdf[4]
print(len(pre_ecdf[2]))
print(pre_ecdf[3])
wins = np.arange(len(pre_ecdf[2]))
wins = np.delete(wins, [int(len(pre_ecdf[2])/2) - 1, int(len(pre_ecdf[2])/2)])


print(obs)

#n_list = [n.neut_ks(arr[x], pre_ecdf[2][x], lins, nn, 1e-6) for x in range(len(pre_ecdf[3]))]
#n_dict = {str(i): n_list[i] for i in range(len(n_list))}
def sym(dict1):
    c = dict1["c"]
    print(f'sbs: {c}')
    list1 = [n.sbs_kl_solve(arr[x], pre_ecdf[2][x], vec1, lins, nn, 0.5, 1, pre_ecdf[3][x], mut, c=c) for x in wins]
    #print(list1)
    return {str(i): list1[i] for i in range(len(list1))}
def over(dict2):
    c = dict2["c"]
    print(f'overdom: {c}')
    list1 = [n.overdom_kl_solve(arr[x], pre_ecdf[2][x], vec1, lins, c, p, pre_ecdf[3][x], mut) for x in wins]
    #print(list1)
    return {str(i): list1[i] for i in range(len(list1))} 
'''
def neut(dict3):
    ne = dict3["ne"]
    n_list = [n.neut_kl(arr[x], np.array(pre_ecdf[5]), pre_ecdf[2][x], lins, ne*nn, mut) for x in wins]
    #print(f'neut list = {n_list}')
    return {str(i): n_list[i] for i in range(len(n_list))}
'''
def dist(dict1, dict2):
    a = np.array(list(dict1.values()))
    b = np.array(list(dict2.values()))
    return np.linalg.norm(a-b)
models = [sym, over]
priors = [pyabc.Distribution(c=pyabc.RV("uniform", 0.6*nn, 0.9*nn)), pyabc.Distribution(c=pyabc.RV("uniform", 0.1*nn, 1.1*nn))]
print(len(priors))
abc = pyabc.ABCSMC(models, priors, pyabc.distance.PNormDistance(p=1))
db_path = "sqlite:///" + os.path.join(tempfile.gettempdir(), "test_6" + db_num + ".db")
bflist=[]
for i in range(3):
    while True:
        try:
            history = abc.new(db_path, obs)
            history = abc.run(minimum_epsilon=0.05, max_nr_populations=12)
            break
        except AssertionError:
            pass

    prob=(history.get_model_probabilities())
    print(prob.iloc[-1])
    bflist.append([prob.iloc[-1][0], prob.iloc[-1][1]])

print(bflist)
arr2 = np.ravel(np.array(bflist))
np.save("results_" + db_num + ".npy", arr2)
