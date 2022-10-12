# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 16:21:43 2022

@author: jihon
"""


import json
import requests
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
from rdkit import Chem
from rdkit.Chem import AllChem


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
    k = np.array([('+' not in s) and ('-' not in s) and ('.' not in s) for s in res['CanonicalSMILES']])
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
        else:
            pass
            # smi = Chem.MolToSmiles(mol, isomericSmiles=False)
        if ('.' in smi) or ('+' in smi) or ('-' in smi):
            continue
        key = res.loc[i, 'InChIKey'].split('-')[0]
        if key not in keys:
            keep.append(i)
            keys.append(key)
    keep, keys = np.array(keep), np.array(keys)
    return res.loc[keep]


def retrieve_by_exact_mass(mass, ppm = 10):
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


def retrieve_by_formula(formula, timeout=999):
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



def retrieve_by_exact_mass_database(mass, database, ppm = 10, priority=[]):
    min_mass = mass - mass * ppm / 10 ** 6
    max_mass = mass + mass * ppm / 10 ** 6
    result = database[np.logical_and(database['Exact mass']>=min_mass, database['Exact mass']<=max_mass)]
    result = result.reset_index(drop=True)
    if len(priority) > 0:
        k = set()
        for pr in priority:
            k = k | set(np.where(result[pr].values.astype(str) != 'nan')[0])
        if len(k) > 0:
            result = result.loc[list(k),:]
    result = result[['Title', 'Formula', 'SMILES', 'InChIkey']]
    if len(result) == 0:
        return []
    result.columns = ['Title', 'MolecularFormula', 'CanonicalSMILES', 'InChIKey']
    result = refine_compound_list(result)
    result = result.reset_index(drop=True)
    return result


def retrieve_by_formula_database(formula, database, priority=[]):   
    result = database[database['Formula'] == formula]
    result = result.reset_index(drop=True)
    if len(priority) > 0:
        k = set()
        for pr in priority:
            k = k | set(np.where(result[pr].values.astype(str) != 'nan')[0])
        if len(k) > 0:
            result = result.loc[list(k),:]
    result = result[['Title', 'Formula', 'SMILES', 'InChIkey']]
    if len(result) == 0:
        return []
    result.columns = ['Title', 'MolecularFormula', 'CanonicalSMILES', 'InChIKey']
    result = refine_compound_list(result)
    result = result.reset_index(drop=True)
    return result



if __name__ == '__main__':
    
    pass