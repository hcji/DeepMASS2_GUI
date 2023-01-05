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
sirius_index = [int(i.split('_')[-1]) for i in sirius_files]

deepmass_path = "example/Plasma/result"
deepmass_files = [name for name in os.listdir(deepmass_path)]
deepmass_index = [int(i.split('_')[-1].split('.')[-2]) for i in deepmass_files]

msfinder_path = "example/Plasma/msfinder/Structure result-2090.txt"
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

    ranking_result.append([name, true_key, sirius_rank, msfinder_rank, deepmass_rank])

ranking_result = pd.DataFrame(ranking_result, columns = ['Challenge', 'True Inchikey2D', 'SIRIUS Ranking', 'MSFinder Ranking',
                                                         'DeepMASS Ranking'])

import seaborn as sns
import matplotlib.pyplot as plt

ratios = []
for i in range(1, 11):
    deepmass_ratio = len(np.where(ranking_result['DeepMASS Ranking'] <= i )[0]) / len(np.where(~np.isnan(ranking_result['DeepMASS Ranking']))[0])
    msfinder_ratio = len(np.where(ranking_result['MSFinder Ranking'] <= i )[0]) / len(np.where(~np.isnan(ranking_result['MSFinder Ranking']))[0])
    sirius_ratio = len(np.where(ranking_result['SIRIUS Ranking'] <= i )[0]) / len(np.where(~np.isnan(ranking_result['SIRIUS Ranking']))[0])
    ratios.append([deepmass_ratio, sirius_ratio, msfinder_ratio])
ratios = pd.DataFrame(ratios, columns = ['DeepMASS', 'SIRIUS', 'MSFinder'])

x = np.arange(1,11)
plt.figure(dpi = 300)
plt.plot(x, ratios['DeepMASS'], label = 'DeepMASS', marker='D', color = '#EE00007F')
plt.plot(x, ratios['SIRIUS'], label = 'SIRIUS', marker='D', color = '#008B457F')
plt.plot(x, ratios['MSFinder'], label = 'MSFinder', marker='D', color = '#3B49927F')
plt.xlim(0.5, 10.5)
plt.ylim(0.4, 0.9)
plt.xticks(np.arange(1, 11, 1))
plt.xlabel('topK', fontsize = 12)
plt.ylabel('ratio', fontsize = 12)
plt.legend(loc='lower right')


'''
import xml.etree.ElementTree as ET
hmdb_xml = ET.parse("E:/HMDB/serum_metabolites.xml")
root = hmdb_xml.getroot()
output = []
for m in tqdm(root):
    [hmdb_id, name, inchikey] = [None, None, None]
    for item in m:
        if item.tag == '{http://www.hmdb.ca}accession':
            hmdb_id = item.text
        if item.tag == '{http://www.hmdb.ca}name':
            name = item.text
        if item.tag == '{http://www.hmdb.ca}inchikey':
            inchikey = item.text
    output.append([hmdb_id, name, inchikey])
output = pd.DataFrame(output, columns = ['HMDB_ID', 'Name', 'Inchikey'])
output.to_csv("E:/HMDB/serum_metabolites.csv", index = False)
'''

hmdb_checking = []
hmdb_compounds = pd.read_csv('E:/HMDB/serum_metabolites.csv')
hmdb_keys = [k[:14] for k in hmdb_compounds['Inchikey'] if type(k) is str]
for s in tqdm(spectrums):
    name = s.metadata['compound_name']
    index = int(name.split('_')[-1])
    
    # rank of SIRIUS
    try:
        sirius_file = "/{}/structure_candidates.tsv".format(sirius_files[sirius_index.index(index)])
        sirius_file = sirius_path + sirius_file
        sirius_result = pd.read_csv(sirius_file, sep='\t')
    except:
        sirius_result = None

    if sirius_result is not None:
        sirius_key = np.array([k for k in sirius_result['InChIkey2D']])
        if len(sirius_key) > 0:
            sirius_top = sirius_key[0]
        else:
            sirius_top = None
    else:
        sirius_top = None
    
    # rank of deepmass
    deepmass_file = "/{}".format(deepmass_files[deepmass_index.index(index)])
    deepmass_file = deepmass_path + deepmass_file
    deepmass_result = pd.read_csv(deepmass_file)
    deepmass_key = np.array([k[:14] for k in deepmass_result['InChIKey']])
    if len(deepmass_key) > 0:
        deepmass_top = deepmass_key[0]
    else:
        deepmass_top = None
        
    # rank of ms-finder
    msfinder_title = np.array([v.replace('unknown', 'unknown_') for v in msfinder_result['Title'].values])
    msfinder_index = np.where(msfinder_title == name)[0]
    if len(msfinder_index) > 0:
        msfinder_key = [str(s)[:14] for s in msfinder_result.loc[msfinder_index[0], msfinder_columns].values]
        msfinder_top = msfinder_key[0]
    else:
        msfinder_top = None

    if sirius_top in hmdb_keys:
        sirius_in = True
    else:
        sirius_in = False
    
    if deepmass_top in hmdb_keys:
        deepmass_in = True
    else:
        deepmass_in = False

        
    if msfinder_top in hmdb_keys:
        msfinder_in = True
    else:
        msfinder_in = False
    hmdb_checking.append([sirius_in, msfinder_in, deepmass_in, ])
    
hmdb_checking = pd.DataFrame(hmdb_checking, columns = ['SIRIUS', 'MSFinder', 'DeepMASS'])

