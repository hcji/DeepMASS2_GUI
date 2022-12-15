# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 09:59:18 2022

@author: DELL
"""

import os
import numpy as np
import pickle
import random

from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem

from matchms.importing import load_from_mgf
from matchms.exporting import save_as_mgf, save_as_msp


ADDUCT_MASSES = {
    # POSITIVE
    '[M+H]+': 1.007276,
    '[M-H2O+H]+': 1.007276 - 18.010565, '[M+H-H2O]+': 1.007276 - 18.010565,
    '[M-2H2O+H]+': 1.007276 - 2 * 18.010565, '[M+H-2H2O]+': 1.007276 - 2 * 18.010565,
    '[M+Na]+': 22.98922,
    '[M]+': 0,
    '[M+NH4]+': 18.03383,
    '[M+H-NH3]+': 1.007276 - 17.026549,
    # NEGATIVE
    '[M-H]-': -1.007276,
    '[M+Cl]-': 34.969402,
    '[M+FA-H]-': 44.9982, '[M-H+FA]-': 44.9982,
    '[M]-': 0,
    '[M-2H]-': -1.007276
}


def get_true_precursor_mz_from_mass(mass, precursor_type):
    """
    Calculate precursor mz based on (exact or monoisotopic) mass and precursor type
    :param mass: scalar, mass of a compound, e.g. monoisotopic or exact mass.
    :param precursor_type: string, precursor type, e.g. '[M+H]+'
    :return: scalar, ion mass / precursor mz
    """
    try:
        return mass + ADDUCT_MASSES[precursor_type]
    except KeyError:
        raise KeyError("Unsupported precursor-type '%s'." % precursor_type)


outfile = os.path.join('example/Test/test_spectrums_positive.pickle')
with open(outfile, 'rb') as file:
    testing_spectrums_ = pickle.load(file)

outfile = os.path.join('example/Test/test_spectrums_negative.pickle')
with open(outfile, 'rb') as file:
    testing_spectrums_ += pickle.load(file)

testing_spectrums = []
for i,s in enumerate(tqdm(testing_spectrums_)):
    mol = Chem.MolFromSmiles(s.metadata['smiles'])
    wt = AllChem.CalcExactMolWt(mol)
    smi = Chem.MolToSmiles(mol, isomericSmiles=False)
    try:
        precursor_mz = get_true_precursor_mz_from_mass(wt, s.get('adduct'))
    except:
        continue
    s = s.set('compound_name', 'challenge_{}'.format(len(testing_spectrums)))
    s = s.set('smiles', smi)
    s = s.set('parent_mass', wt)
    s = s.set('precursor_mz', precursor_mz)
    testing_spectrums.append(s)

# save as mgf
save_as_mgf(list(testing_spectrums), 'example/Test/all_test_set.mgf')

# save as msp individually
testing_spectrums = [s for s in load_from_mgf('example/Test/all_test_set.mgf')]
for i, s in enumerate(testing_spectrums):
    path = 'example/Test/msp/challenge_{}.msp'.format(i)
    save_as_msp([testing_spectrums[i]], path)
    
    with open(path, encoding = 'utf-8') as msp:
        lines = msp.readlines()
        lines = [l.replace('_', '') for l in lines]
        lines = [l.replace('ADDUCT', 'PRECURSORTYPE') for l in lines]
    with open(path, 'w') as msp:
        msp.writelines(lines)

    # only for ms-finder
    path_msfinder = path.replace('/msp/', '/msfinder/')
        
    # exclude large compound, as ms-finder is very slow for them
    if float(s.metadata['parent_mass']) >= 850:
        continue
    with open(path_msfinder, 'w') as msp:
        msp.writelines(lines)