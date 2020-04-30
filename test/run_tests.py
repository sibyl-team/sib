#!/usr/bin/env python3
## Author: Fabio Mazza
import unittest
import sys
import numpy as np
import pandas as pd
import os
import data_load

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(1,dir_path+"/..")
import sib
from pathlib import Path

script_path = Path(dir_path)

FOLDER = script_path / "data_tree"
NUM_CPUS=10

def callback(t,err,f):
    print(f"{t:4d}, {err:3.2e}",end="\r")

class SibillaTest(unittest.TestCase):

    def find_sources_sib(self,obs_sib,full_epi):
        mu = self.params["mu"]
        src = np.where(full_epi[0])[0][0]
        sib_pars = sib.Params(prob_r=sib.Gamma(mu=mu))
        sib_fg = sib.FactorGraph(sib_pars,self.contacts_sib,obs_sib)

        sib.iterate(sib_fg,maxit=1000,tol=6e-6,callback=callback)
        
        #print("\n",end="")
        #iterate_damp(sib_fg,2000,callback,0.5)

        p_sources = []
        for n in sib_fg.nodes:
            marg = sib.marginal(n)
            p_sources.append(marg[1])

        p_sources = np.array(p_sources)
        nodes = np.argsort(p_sources[:,1])[::-1]
        accu = np.cumsum(nodes==src)
        return p_sources,accu

    def setUp(self):
        # params,contacts,observ,all_epi
        self.data= data_load.load_exported_data(FOLDER)
        self.params = self.data[0]
        df_cont = self.data[1][["i","j","t","lambda"]]
        self.contacts_sib=list(zip(*[df_cont[k] for k in df_cont.keys()]))
        self.obs_all_sib = []
        for ob in self.data[2]:
            df_obs = data_load.convert_obs_to_df(ob)
            df_obs = df_obs[["i","st","t"]]
            obs_sib = list(df_obs.to_records(index=False))
            self.obs_all_sib.append(obs_sib)

        sib.set_num_threads(NUM_CPUS)
    
    def test_inference(self):
        print("Executing run 1")
        probs1 = np.stack([self.find_sources_sib(obs,epi)[0] for obs,epi in zip(self.obs_all_sib,self.data[3])])
        #print(probs1[3][0])
        self.assertEqual(np.any(probs1 == np.inf),False)

        print("\nExecuting run 2")
        probs2 = np.stack([self.find_sources_sib(obs,epi)[0] for obs,epi in zip(self.obs_all_sib,self.data[3])])

        print("")
        self.assertEqual(np.all( (probs1-probs2) < 1e-7 ),True)

    def test_accuracy(self):
        print("Test accuracy")
        accu_all = np.stack([self.find_sources_sib(obs,epi)[1] for obs,epi in zip(self.obs_all_sib,self.data[3])])
        accu_curve = accu_all.mean(0)

        accu_meas = accu_curve.cumsum().sum()/self.params["n"]
        self.assertGreaterEqual(accu_meas,12,msg="The accuracy does not correspond. Maybe the observations have the wrong order?")
    
    #def test_

if __name__ == '__main__':
    unittest.main(failfast=True)