# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 15:27:35 2022

@author: DELL
"""

import os
import numpy as np
import pandas as pd
import pubchempy as pcp

from tqdm import tqdm
from matchms import Spectrum
from matchms.importing import load_from_mgf
from matchms.exporting import save_as_mgf, save_as_msp
from core.identification import spectrum_processing
from core.msdial import load_MS_DIAL_Alginment

ident_pos = pd.read_excel("E:/Data/tomato_data/MH2102-009-DR-TMT-pos-LS-FINAL-MSMS.xlsx").loc[:,['Compound Name', 'm/z (Expected)', 'RT']]
ident_neg = pd.read_excel("E:/Data/tomato_data/MH2102-009-TMT-neg-LS.xlsx").loc[:,['Compound Name', 'm/z (Expected)', 'RT']]

ident_pos['Ionization'] = 'POS'
ident_neg['Ionization'] = 'NEG'

def remove_uncertain(data):
    keep = []
    for i in data.index:
        mz = data.loc[i, 'm/z (Expected)']
        rt = data.loc[i, 'RT']
        if np.isnan(mz) or np.isnan(rt):
            continue
        
        check = np.logical_and(np.abs(data.loc[:, 'm/z (Expected)'].values - mz) < 0.05, 
                               np.abs(data.loc[:, 'RT'].values - rt) < 0.5)
        if len(np.where(check)[0]) > 1:
            continue
        else:
            keep.append(i)
            
    final = data.loc[keep,:]
    final = final.reset_index(drop = True)
    return final

ident_pos = remove_uncertain(ident_pos)
ident_neg = remove_uncertain(ident_neg)
tomato_identified = ident_pos.append(ident_neg, ignore_index = True)

'''
exact_mass, inchikey, molecular_formula = [], [], []
for i in tqdm(tomato_identified.index):
    name = tomato_identified.loc[i, 'Compound Name']
    if name[-1] == ' ':
        name = name[:-1]
    compounds = pcp.get_compounds(name, 'name')

    m = np.nan
    k = np.nan
    f = np.nan
    if len(compounds) > 0:
        cmd = compounds[0]
        m = float(cmd.exact_mass)
        k = cmd.inchikey
        f = cmd.molecular_formula
    else:
        name = name.replace('-', '')
        compounds = pcp.get_compounds(name, 'name')
        if len(compounds) > 0:
            cmd = compounds[0]
            m = float(cmd.exact_mass)
            k = cmd.inchikey
            f = cmd.molecular_formula
    
    exact_mass.append(m)
    inchikey.append(k)
    molecular_formula.append(f)

tomato_identified['exact_mass'] = exact_mass
tomato_identified['inchikey'] = inchikey
tomato_identified['molecular_formula'] = molecular_formula
tomato_identified.to_csv('example/Tomato/tomato_identified.csv', index = False)
'''

