"""
Library for saving and loading instances in the benchmark repository

Author: Fabio Mazza


"""
import numpy as np
import pandas as pd

import json
import bz2
import tarfile
from pathlib import Path

def save_json(obj,file_,indent=1):
    f = open(file_,"w")
    json.dump(obj,f,indent=indent)
    f.close()


CONTACTS_FNAME = "contacts.csv.bz2"
EPI_NAME = "epidemy"
OBS_NAME = "observ"
OBS_JSON_FNAME = "obs.json"
PARAMS_FNAME = "parameters.json"
ALL_OBS_FILE= "all_obs.json.bz2"
ALL_EPI = "all_epidemies"

NAMES_COLS_CONTACTS = ["t","i","j","lambda"]
CONTACTS_DTYPES = dict(zip(NAMES_COLS_CONTACTS,(np.int,np.int,np.int,np.float) ))



def load_exported_data(folder_path,epidemies_with_name=False,pandas_df=True):
    """
    Load only the binary formatted files
    """
    fold = Path(folder_path)
    if not fold.exists():
        raise ValueError("Folder doesn't exists")

    with open(fold / PARAMS_FNAME) as f:
        params = json.load(f)

    with bz2.open(fold/(CONTACTS_FNAME),"r") as f:
        if pandas_df:
            contacts = pd.read_csv(f,names=NAMES_COLS_CONTACTS,dtype=CONTACTS_DTYPES)
        else:
            contacts = np.loadtxt(f,delimiter=",")

    file_all_obs = fold / ALL_OBS_FILE
    if file_all_obs.exists():
        with bz2.open(file_all_obs,"rt") as f:
            observ = json.load(f)
    else:
        print("No observations found")
        observ = None
    
    all_epi = np.load(fold/(ALL_EPI+".npz"))
    epids = all_epi.files
    indc  = [int(n.split("_")[1]) for n in epids]
    r = list(zip(indc,epids))
    r.sort()    
    epid_stacked = np.stack([all_epi[name] for i,name in r])

    return params,contacts,observ,epid_stacked

def convert_obs_to_df(globs):
    import itertools
    STATE_MAP = {"S":0,"I":1,"R":2}
    columns=["st","i","t"]
    
    all_dfs = []
    for state,obse in globs.items():
        v = STATE_MAP[state]
        for t, nodes in obse.items():
            data_iter = itertools.product([v],nodes,[int(t)])
            all_dfs.append(pd.DataFrame(data_iter,columns=columns))
    return pd.concat(all_dfs)
    