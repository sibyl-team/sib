#!/usr/bin/env python3
# Copyright 2020 Sybil team
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
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
FIELDS_FILE = "fields_tree.npz"


def callback(t, err, f):
    print(f"{t:4d}, {err:3.2e}", end="\r")




def load_run_data(filename,n_inst,num_nodes):
    data = np.load(filename)
    all_data = []
    for inst in range(n_inst):
        margist = [data[f"{inst}_{n}"] for n in range(num_nodes)]
        all_data.append(margist)
    return all_data




class SibillaTest(unittest.TestCase):

    def run_sib_instance(self,inst,callback_fun=callback):
        mu = self.params["mu"]
        sib_pars = sib.Params(prob_r=sib.Gamma(mu=mu))
        sib_fg = sib.FactorGraph(sib_pars, self.contacts_sib, self.obs_all_sib[inst])

        sib.iterate(sib_fg, maxit=1000, tol=6e-6, callback=callback_fun)

        return sib_fg

    def calc_beliefs_fields(self,inst):
        fg = self.run_sib_instance(inst)

        data = []
        for n in fg.nodes:
            res = (np.array(n.bt),np.array(n.bg))
            fields = (np.array(n.ht),np.array(n.hg))
            data.append((np.stack(res),np.stack(fields)))
        return data

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
        self.loaded_beliefs = load_run_data(script_path/BELIEFS_FILE,self.n_inst,self.num_nodes)
        self.loaded_fields = load_run_data(script_path/FIELDS_FILE,self.n_inst,self.num_nodes)

    def test_inference(self):
        print("\n--- Executing trial runs ---")
        print("Run 1")
        probs1 = np.stack([self.find_sources_sib(i)[0] for i in range(self.n_inst)])
        # print(probs1[3][0])
        self.assertEqual(np.any(probs1 == np.inf), False)

        print("\nRun 2")
        probs2 = np.stack([self.find_sources_sib(i)[0] for i in range(self.n_inst)])

        print("")
        self.assertEqual(np.all((probs1-probs2) < 1e-12), True)

    def test_accuracy(self):
        print("\n--- Testing accuracy ---")
        accu_all = np.stack([self.find_sources_sib(i)[1] for i in range(self.n_inst) ])
        accu_curve = accu_all.mean(0)

        accu_meas = accu_curve.cumsum().sum()/self.params["n"]
        self.assertGreaterEqual(
            accu_meas, 12, msg="The accuracy does not correspond. Maybe the observations have the wrong order?")

    def test_beliefs_fields(self):
        print("\n--- Testing beliefs and fields ---")
        beliefs_fields = [self.calc_beliefs_fields(i) for i in range(self.n_inst)]
        
        for i in range(self.n_inst):
            for n in range(self.num_nodes):
                #with self.subTest(inst=i,node=n):
                msg_fail = "Test on inst {} for node {} failed".format(i,n)
                
                self.assertEqual(np.all(beliefs_fields[i][n][0] - self.loaded_beliefs[i][n] < 1e-12),True,msg_fail)
                
                self.assertEqual(np.all(beliefs_fields[i][n][1] == self.loaded_fields[i][n]),True,msg_fail)


if __name__ == '__main__':
    unittest.main(failfast=False)
