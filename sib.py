from _sib import *

def marginal(n):
    rbt = [0.0] * len(n.bt)
    lbg = [0.0] * len(n.bg)
    T = len(n.bt) - 1
    lbg[0] = n.bg[0]
    rbt[T] = n.bt[T]
    for t in range(1, T+1):
        lbg[t] = lbg[t-1] + n.bg[t]
        rbt[T - t] = rbt[T - t + 1] + n.bt[T - t]
    marg=[]
    for t in range(1, T):
        marg.append((rbt[t], 1-rbt[t]-lbg[t-1], lbg[t-1]))
    return marg

def iterate(f, maxit=100, tol=1e-3, damping=0.0, callback=(lambda t, err, f: print(t, err, flush=True))):
    for t in range(maxit):
        err = f.update(damping)
        if callback(t, err, f) == False:
            break;
        if err < tol:
            break;
