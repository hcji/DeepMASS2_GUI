# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 09:51:53 2023

@author: DELL
"""


import json

import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
from rdkit import Chem
from rdkit.Chem import AllChem


def search_from_database(s, database, ppm = 10):
    mass = s.get('parent_mass')
    formula = s.get('formula')
    if formula is not None:
        candidates = retrieve_by_formula_database(formula, database)
    elif mass is not None:
        candidates = retrieve_by_exact_mass_database(mass, database, ppm = ppm)
    else:
        candidates = None
    if candidates is not None:
        candidates['Formula Score'] = np.nan
        candidates['Structure Score'] = np.nan
        candidates['Consensus Score'] = np.nan
    s.set('annotation', candidates)
    return s


def search_from_pubchem(s, ppm = 10):
    mass = s.get('parent_mass')
    formula = s.get('formula')
    if formula is not None:
        try:
            candidates = retrieve_by_formula_pubchem(formula)
        except:
            candidates = None
    elif mass is not None:
        try:
            candidates = retrieve_by_exact_mass_pubchem(mass, ppm = ppm)
        except:
            candidates = None
    else:
        candidates = None
    if candidates is not None:
        candidates['Formula Score'] = np.nan
        candidates['Structure Score'] = np.nan
        candidates['Consensus Score'] = np.nan
    s.set('annotation', candidates)
    return s


def retrieve_by_cid_list(idlist):
    res = []
    idstring = ''
    for i, cid in enumerate(idlist):
        idstring += ',' + str(cid)
        if ((i%100==99) or (i==len(idlist)-1)):
            url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + idstring[1:(len(idstring))] + "/property/MolecularFormula,InChIKey,CanonicalSMILES/JSON"
            soup = BeautifulSoup(requests.get(url, timeout=15).content, "html.parser")
            properties = json.loads(str(soup))['PropertyTable']['Properties']
            properties = [[s[k] for k in s.keys()] for s in properties]
            res += properties
            idstring = ''
    res = pd.DataFrame(res)
    res.columns = ['CID','MolecularFormula','CanonicalSMILES','InChIKey']
    k = np.array([('.' not in s) for s in res['CanonicalSMILES']])
    res = res.loc[k,:]
    res = res.reset_index(drop = True)
    return res


def refine_compound_list(res):
    keep, keys = [], []
    for i in res.index:
        smi = res.loc[i, 'CanonicalSMILES']
        formula = res.loc[i, 'MolecularFormula']
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        formula_cal = AllChem.CalcMolFormula(mol)
        if formula != formula_cal:
            continue
        if '.' in smi:
            continue
        key = res.loc[i, 'InChIKey'].split('-')[0]
        if key not in keys:
            keep.append(i)
            keys.append(key)
    keep, keys = np.array(keep), np.array(keys)
    return res.loc[keep]


def retrieve_by_exact_mass_pubchem(mass, ppm = 10):
    min_mass = mass - mass * ppm / 10 ** 6
    max_mass = mass + mass * ppm / 10 ** 6
    url = '''https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&retmode=json&term={}%3A{}[ExactMass]&retmax=9999'''.format(min_mass, max_mass)
    try:
        res = requests.get(url).text
    except:
        raise ConnectionError
    res = json.loads(res)
    idlist = res['esearchresult']['idlist']
    result = retrieve_by_cid_list(idlist)
    result = refine_compound_list(result)
    result = result.reset_index(drop = True)
    return result


def retrieve_by_formula_pubchem(formula, timeout=999):
    url = '''https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/{}/cids/json'''.format(formula)
    try:
        res = requests.get(url).text
    except:
        raise ConnectionError
    res = json.loads(res)
    idlist = res['IdentifierList']['CID']
    result = retrieve_by_cid_list(idlist)
    result = refine_compound_list(result)
    result = result.reset_index(drop = True)
    return result


def retrieve_by_exact_mass_database(mass, database, ppm = 10):
    min_mass = mass - mass * ppm / 10 ** 6
    max_mass = mass + mass * ppm / 10 ** 6
    result = database[np.logical_and(database['Exact mass']>=min_mass, database['Exact mass']<=max_mass)]
    result = result[['Title', 'Formula', 'SMILES', 'InChIkey', 'Database IDs']]
    if len(result) == 0:
        return None
    result.columns = ['Title', 'MolecularFormula', 'CanonicalSMILES', 'InChIKey', 'Database IDs']
    result = refine_compound_list(result)
    result = result.reset_index(drop=True)
    return result


def retrieve_by_formula_database(formula, database):   
    result = database[database['Formula'] == formula]
    result = result[['Title', 'Formula', 'SMILES', 'InChIkey', 'Database IDs']]
    if len(result) == 0:
        return None
    result.columns = ['Title', 'MolecularFormula', 'CanonicalSMILES', 'InChIKey', 'Database IDs']
    result = refine_compound_list(result)
    result = result.reset_index(drop=True)
    return result



if __name__ == '__main__':
    
    '''
    from core.importing.load_from_files import load_from_files
    
    database = pd.read_csv('data/DeepMassStructureDB-v1.1.csv')
    spectrums = load_from_files(["example/all_casmi.mgf"])
    spectrums = [search_from_database(s, database, ppm = 10) for s in spectrums]
    '''
