#!/usr/bin/env python
# coding: utf-8
import sys
sys.path.insert(0,"..")

import numpy as np
import json
import sib

with open("data/large_tree_pars.json") as f:
    params = json.load(f)

cts_beliefs_f = np.load("data/large_tree_data.npz")

cts = cts_beliefs_f["cts"]

beliefs_all  = cts_beliefs_f["beliefs"]

obs_all = {}
for k in cts_beliefs_f.files:
    if "obs_" in k:
        u=int(k.split("_")[-1])
        #print(k, u)
        obs_all[u] = cts_beliefs_f[k]

#close file
cts_beliefs_f.close()

sib_pars = sib.Params(prob_r=sib.Gamma(mu=params["mu"]))

N = params["N"]
t_limit = params["t_limit"]
cts_sib = [(int(r["i"]),int(r["j"]),int(r["t"]),r["lam"]) for r in cts]

tests = [sib.Test(s==0,s==1,s==2) for s in range(3)]
def make_obs_sib(N, t_limit,obs, tests):
    obs_list_sib =[(i,-1,t) for t in [t_limit] for i in range(N) ]
    obs_list_sib.extend([(r["i"],tests[r["st"]],r["t"]) for r in obs])

    obs_list_sib.sort(key=lambda x: x[-1])

    return obs_list_sib

callback = lambda t, err, fg: print(f"iter: {t:6}, err: {err:.5e} ", end="\r")

for ii,obs in obs_all.items():
    fg = sib.FactorGraph(params=sib_pars)
    beliefs = beliefs_all[ii]
    #print(f"Instance {ii}")
    #for c in cts_sib:
    #    fg.append_contact(*c)
    fg.append_contacts_npy(cts["i"], cts["j"], cts["t"], cts["lam"])
    obs_list_sib = make_obs_sib(N,t_limit, obs, tests)
    for o in obs_list_sib:
        fg.append_observation(*o)
    
    sib.iterate(fg,200,1e-20,callback=callback )
    print("")
    s=0.
    for i in range(len(fg.nodes)):
        s+=np.abs(np.array(fg.nodes[i].bt)-beliefs[i][0]).sum()
        s+=np.abs(np.array(fg.nodes[i].bg)-beliefs[i][1]).sum()
        #fg.nodes[i].bg]))
    print(f"instance {ii}: {s:4.3e} {s < 1e-10}")





