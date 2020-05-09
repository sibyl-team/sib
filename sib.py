from _sib import *

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
        n = f.nodes[f.get_index(i)]
        ttrue = list(n.times).index(t)-1
        M[i] = n.marginal_t(ttrue)
        #sib.marginal

    return M

def marginal_t(i, t):
    '''
    returns the marginals of nodes at fixed time 
    
    - i: name node
    - t: time
    
    return: list - probability to be [prob_S, prob_I, prob_R]
    '''
    n = f.nodes[f.get_index(i)]
    ttrue = list(n.times).index(t)-1
    M = n.marginal_t(ttrue)
    #sib.marginal

    return M



def iterate(f, maxit=100, tol=1e-3, damping=0.0, callback=(lambda t, err, f: print(t, err, flush=True))):
    for t in range(maxit):
        err = f.update(damping)
        if callback(t, err, f) == False:
            break;
        if err < tol:
            break;
