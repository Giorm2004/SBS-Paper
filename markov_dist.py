import numpy as np
from scipy import stats
from scipy import linalg


s1 =np.ones((3,1))


def g_mat(n, p, r):
    c = 1/(2*n*p)
    w = - c - 0.5*r 
    return np.matrix([[w, 0.5*r, 0], [0.5*r,  -1*r, 0.5*r], [0, 0.5*r, w]])
def state_vec(s,g):
    p = (4+3*s)/(4+s)
    p0 = 1/(1+(p**(-1*g/2)))
    return np.array([[p0**2, 2*(1-p0)*p0, (1-p0)**2]])
def ph_cdf(mat, state_vec1, x):
    mat = mat.reshape(3,3)
    return 1 - state_vec1.dot(linalg.expm(x*mat).dot(s1))

def ph_pdf(mat, state_vec1,x):
    mat = mat.reshape(3,3)
    s = -1*mat.dot(s1)
    return state_vec1.dot(linalg.expm(x*mat).dot(s))
def ex_t(mat, state_vec1):
    U =linalg.inv(-1*mat)
    u1= U.dot(s1)
    return state_vec1.dot(u1)[0,0]
mat = g_mat(1e4, 0.5, 1e-4)
vec = state_vec(0.5, 10)
class coal_time_dist(stats.rv_continuous):
    def _argcheck(self, mat, vec):
        return (type(mat) == np.matrix) & (type(vec) == np.ndarray)
    def _cdf(self,x, mat, vec):
        x = float(x[0])
        mat=mat.reshape(3,3)
        vec = vec[0:3,]
        return ph_cdf(mat, vec, x)[0]
    def _pdf(self,x, mat, vec):
        x = float(x[0])
        mat=mat.reshape(3,3)
        vec = vec[0:3,]
        return ph_pdf(mat, vec, x)[0]
    def _rvs(self, mat, vec, size=1, random_state=None):
        umax = 5*ph_pdf(mat, vec,ex_t(mat, vec))[0,0]
        vmax = 5*ex_t(mat, vec)
        def f(x):
            return self.pdf(x, mat, vec)
        gen = stats.sampling.RatioUniforms(f, umax=umax, vmin=0, vmax=vmax, random_state=random_state)
        return gen.rvs(size)



dist = coal_time_dist(name="dist", badvalue=6)

