from _sib import *


def iterate(f, maxit=100, tol=1e-3, callback=(lambda t, err, f: print(t, err))):
    for t in range(maxit):
        err = f.update()
        if callback(t, err, f) == False
            break;
        if err < tol:
            break;
