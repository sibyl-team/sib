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


def FactorGraph(params = _sib.Params(_sib.Uniform(1.0), _sib.Exponential(0.5), 0.1, 0.45, 0.0, 0.0 ,0.0, 0.0),
                contacts = [],
                observations = [],
                tests = [],
                times = [],
                individuals = []):
    if len(observations) > 0:
        if len(tests) > 0:
            print("only one between tests and observations is allowed")
            return None
        tests = [(i, params.fakeobs if s == -1 else params.obs[s],t) for (i,s,t) in observations if s <= 2]
    return _sib.FactorGraph(params = params, contacts = contacts, tests = tests, individuals = individuals)


def iterate(f,
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
        err = f.update(damping=damping, learn=learn)
        if callback(t, err, f)  == False:
            break;
        if err < tol:
            break;
    if newline:
        print()
