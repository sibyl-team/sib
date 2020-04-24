from _sib import *


def iterate(f, maxit=100, tol=1e-3, callback=(lambda t, err, f: print(t, err))):
    for t in range(maxit):
        err = f.update()
        callback(t, err, f)
        if err < tol:
            break;
