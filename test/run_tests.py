#!/usr/bin/env python3
# This file is part of sibilla : inference in epidemics with Belief Propagation
# Author: Fabio Mazza

from pathlib import Path
import unittest
import sys
import numpy as np
import pandas as pd
import os
import data_load

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(1, dir_path+"/..")
import sib
script_path = Path(dir_path)

FOLDER = script_path / "data_tree"
NUM_CPUS = 10

BELIEFS_FILE = "beliefs_tree.npz"


def callback(t, err, f):
    print(f"{t:4d}, {err:3.2e}", end="\r")




def load_beliefs(filename,n_inst,num_nodes):
    data = np.load(filename)
    all_marg_load = []
    for inst in range(n_inst):
        margist = [data[f"{inst}_{n}"] for n in range(num_nodes)]
        all_marg_load.append(margist)
    return all_marg_load


class SibillaTest(unittest.TestCase):

    def run_sib_instance(self,inst,callback_fun=callback):
        mu = self.params["mu"]
        sib_pars = sib.Params(prob_r=sib.Gamma(mu=mu))
        sib_fg = sib.FactorGraph(sib_pars, self.contacts_sib, self.obs_all_sib[inst])

        sib.iterate(sib_fg, maxit=1000, tol=6e-6, callback=callback_fun)

        return sib_fg

    def calc_beliefs(self,inst):
        fg = self.run_sib_instance(inst)

        beliefs = []
        for n in fg.nodes:
            res = (np.array(n.bt),np.array(n.bg))
            beliefs.append(np.stack(res))
        return beliefs

    def find_sources_sib(self,inst):
        src = self.sources[inst]
        sib_fg = self.run_sib_instance(inst)

        p_sources = []
        for n in sib_fg.nodes:
            marg = sib.marginal(n)
            p_sources.append(marg[1])

        p_sources = np.array(p_sources)
        nodes = np.argsort(p_sources[:, 1])[::-1]
        accu = np.cumsum(nodes == src)
        return p_sources, accu

    def setUp(self):
        # params,contacts,observ,all_epi
        self.data = data_load.load_exported_data(FOLDER)
        self.params = self.data[0]
        # Format contacts
        df_cont = self.data[1][["i", "j", "t", "lambda"]]
        self.contacts_sib = list(zip(*[df_cont[k] for k in df_cont.keys()]))
        # Format observations
        self.obs_all_sib = []
        for ob in self.data[2]:
            df_obs = data_load.convert_obs_to_df(ob)
            df_obs = df_obs[["i", "st", "t"]]
            obs_sib = list(df_obs.to_records(index=False))
            self.obs_all_sib.append(obs_sib)
        # Find sources nodes
        self.sources = []
        for epi in self.data[3]:
            src = np.where(epi[0])[0][0]
            self.sources.append(src)
        
        self.n_inst = len(self.obs_all_sib)
        self.num_nodes = self.params["n"]

        sib.set_num_threads(NUM_CPUS)
        ## LOAD BELIEFS
        self.loaded_beliefs = load_beliefs(script_path/BELIEFS_FILE,self.n_inst,self.num_nodes)

    def test_inference(self):
        print("\n--- Executing trial runs ---")
        print("Run 1")
        probs1 = np.stack([self.find_sources_sib(i)[0] for i in range(self.n_inst)])
        # print(probs1[3][0])
        self.assertEqual(np.any(probs1 == np.inf), False)

        print("\nRun 2")
        probs2 = np.stack([self.find_sources_sib(i)[0] for i in range(self.n_inst)])

        print("")
        self.assertEqual(np.all((probs1-probs2) < 1e-7), True)

    def test_accuracy(self):
        print("\n--- Testing accuracy ---")
        accu_all = np.stack([self.find_sources_sib(i)[1] for i in range(self.n_inst) ])
        accu_curve = accu_all.mean(0)

        accu_meas = accu_curve.cumsum().sum()/self.params["n"]
        self.assertGreaterEqual(
            accu_meas, 12, msg="The accuracy does not correspond. Maybe the observations have the wrong order?")

    def test_beliefs(self):
        print("\n--- Testing beliefs ---")
        new_all_beliefs = [self.calc_beliefs(i) for i in range(self.n_inst)]
        for i in range(self.n_inst):
            for n in range(self.num_nodes):
                self.assertEqual(np.all(new_all_beliefs[i][n] == self.loaded_beliefs[i][n]),True)


if __name__ == '__main__':
    unittest.main(failfast=True)
