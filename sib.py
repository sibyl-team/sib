from _sib import *

def marginal(n):
    return n.marginal()

def iterate(f, maxit=100, tol=1e-3, damping=0.0, callback=(lambda t, err, f: print(t, err, flush=True))):
    for t in range(maxit):
        err = f.update(damping)
        if callback(t, err, f) == False:
            break;
        if err < tol:
            break;
