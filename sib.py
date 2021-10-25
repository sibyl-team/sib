from _sib import *

import sys

import _sib

__module_file__ = _sib.__file__

def marginal(n):
    return n.marginal()

def marginals_t(f, t):
    '''
    returns the marginals of nodes at fixed time 
    
    - f: sib.f class function
    - t: time
    
    return: dict - probability to be {i : [prob_S, prob_I, prob_R]
    '''
    M = {}
    for i in range(len(f.nodes)):
        M[i] = marginal_t(f.nodes[i], t)

    return M

def marginal_t(n, t):
    '''
    returns the marginals of nodes at fixed time 
    
    - n: sib.Node class
    - t: time
    
    return: list - probability to be [prob_S, prob_I, prob_R]
    '''
    
    # we use "-1", marginal_index removes the source times.
    ttrue = list(n.times).index(t)-1
    M = n.marginal_index(ttrue)

    return M




class FactorGraph(_sib.FactorGraph):


    def gettest(self, s):
        if isinstance(s, int) and s == -1:
            return self._fakeobs
        elif isinstance(s, int) and 0 <= s < len(self.puretest):
            return self.puretest[s]
        else:
            return s

    def __init__(self, params = _sib.Params(_sib.Uniform(1.0), _sib.Exponential(0.5), 0.1, 0.45, 0.0, 0.0),
                contacts = [],
                observations = [],
                tests = [],
                times = [],
                individuals = []):
        self.puretest = [_sib.Test(s==0,s==1,s==2) for s in range(3)]
        self._fakeobs = _sib.Test(1,1,1)
        tests = [(i, self.gettest(s), t) for (i,s,t) in observations+tests]
        _sib.FactorGraph.__init__(self, params = params, contacts = contacts, tests = tests, individuals = individuals)

    def append_observation(self, i, s, t):
        _sib.FactorGraph.append_observation(self, i, self.gettest(s), t)

    def iterate(self,
            maxit=100,
            tol=1e-3,
            damping=0.0,
            learn=False,
            callback=False
        ):
        newline = False
        if callback == False:
            callback = lambda t,e,f : print(f"sib.iterate(damp={damping}): {t}/{maxit} {e:1.3e}/{tol}", end='      \r', flush=True)
            newline = True
        if callback == None:
            callback = lambda t,e,f : None
        for t in range(maxit):
            err = self.update(damping=damping, learn=learn)
            if callback(t, err, self)  == False:
                break;
            if err < tol:
                break;
        if newline:
            print()

def iterate(f, maxit=100, tol=1e-3, damping=0.0, learn=False, callback=False):
    return f.iterate(maxit=maxit, tol=tol, damping=damping, learn=learn, callback=callback)

