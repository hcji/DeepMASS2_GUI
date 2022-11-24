# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 09:13:11 2022

@author: DELL
"""

import os
import numpy as np
import pandas as pd

from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem, inchi
from matchms.importing import load_from_mgf

spectrums = [s for s in load_from_mgf('example/Plasma/ms_ms_plasma.mgf')]

sirius_path = "example/Plasma/sirius"
sirius_files = [name for name in os.listdir(sirius_path) if os.path.isdir(os.path.join(sirius_path, name)) ]
sirius_index = [int(i.split('_')[-2]) for i in sirius_files]

deepmass_path = "example/Plasma/result_all"
deepmass_files = [name for name in os.listdir(deepmass_path)]
deepmass_index = [int(i.split('_')[-1].split('.')[-2]) for i in deepmass_files]

deepmass_path2 = "example/Plasma/result_silicon"
deepmass_files2 = [name for name in os.listdir(deepmass_path2)]
deepmass_index2 = [int(i.split('_')[-1].split('.')[-2]) for i in deepmass_files2]

msfinder_path = "example/Plasma/msfinder/Structure result-2086.txt"
msfinder_result = pd.read_csv(msfinder_path, sep = '\t')
msfinder_columns = [col for col in msfinder_result.columns if 'InChIKey' in col]

ranking_result = []
for s in tqdm(spectrums):
    name = s.metadata['compound_name']
    index = int(name.split('_')[-1])
    true_key = s.get('inchikey')
    if true_key is None:
        continue
    else:
        true_key = true_key[:14]
    
    # rank of SIRIUS
    try:
        sirius_file = "/{}/structure_candidates.tsv".format(sirius_files[sirius_index.index(index)])
        sirius_file = sirius_path + sirius_file
        sirius_result = pd.read_csv(sirius_file, sep='\t')
    except:
        sirius_result = None

    # rank of sirius
    if sirius_result is not None:
        sirius_n = len(sirius_result)
        sirius_key = np.array([k for k in sirius_result['InChIkey2D']])
        sirius_rank = np.where(sirius_key == true_key)[0]
        if len(sirius_rank) == 0:
            sirius_rank = float('inf')
        else:
            sirius_rank = sirius_rank[0] + 1
    else:
        sirius_n = 0
        sirius_rank = float('inf')
    
    # rank of deepmass
    deepmass_file = "/{}".format(deepmass_files[deepmass_index.index(index)])
    deepmass_file = deepmass_path + deepmass_file
    deepmass_result = pd.read_csv(deepmass_file)
    deepmass_key = np.array([k[:14] for k in deepmass_result['InChIKey']])
    deepmass_n = len(deepmass_key)
    deepmass_rank = np.where(deepmass_key == true_key)[0]
    if len(deepmass_rank) == 0:
        deepmass_rank = float('inf')
    else:
        deepmass_rank = deepmass_rank[0] + 1
        
    # rank of deepmass in silicon
    deepmass_file2 = "/{}".format(deepmass_files2[deepmass_index2.index(index)])
    deepmass_file2 = deepmass_path2 + deepmass_file2
    deepmass_result2 = pd.read_csv(deepmass_file2)
    deepmass_key2 = np.array([k[:14] for k in deepmass_result2['InChIKey']])
    deepmass_n2 = len(deepmass_key2)
    deepmass_rank2 = np.where(deepmass_key2 == true_key)[0]
    if len(deepmass_rank2) == 0:
        deepmass_rank2 = float('inf')
    else:
        deepmass_rank2 = deepmass_rank2[0] + 1
        
    # rank of ms-finder
    msfinder_title = np.array([v.replace('unknown', 'unknown_') for v in msfinder_result['Title'].values])
    msfinder_index = np.where(msfinder_title == name)[0]
    if len(msfinder_index) > 0:
        msfinder_key = [str(s)[:14] for s in msfinder_result.loc[msfinder_index[0], msfinder_columns].values]
        msfinder_rank = np.where(np.array(msfinder_key) == true_key)[0]
        if len(msfinder_rank) == 0:
            msfinder_rank = float('inf')
        else:
            msfinder_rank = msfinder_rank[0] + 1
    else:
        msfinder_rank = np.nan

    ranking_result.append([name, true_key, sirius_rank, msfinder_rank, deepmass_rank, deepmass_rank2])

ranking_result = pd.DataFrame(ranking_result, columns = ['Challenge', 'True Inchikey2D', 'SIRIUS Ranking', 'MSFinder Ranking',
                                                         'DeepMASS All Ranking', 'DeepMASS InSilicon Ranking'])

