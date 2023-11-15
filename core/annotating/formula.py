# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 10:16:57 2023

@author: DELL
"""


import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs, AllChem

from core.pycdk import IsotopeFromString, IsotopeSimilarity


def calc_isotopic_score(s):
    if s.get('candidates') is None:
        return s
    elif s.get('isotopic_pattern') is None:
        return s
    else:
        pass
    
    

if __name__ == '__main__':

    from core.Spectrum import Spectrum
    from core.annotating.candidates import search_from_database
    
    s = Spectrum(mz = np.array([]),
                 intensities = np.array([]),
                 isotopic_mz = np.array([222.0230, 223.0264, 224.0188, 225.0220]),
                 isotopic_intensities = np.array([14858164.0, 1567580.875, 666319.625, 71849.90625]),
                 metadata = {'parent_mass': 223.0303,
                             'precursor_mz': 222.0230})
    database = pd.read_csv('data/DeepMassStructureDB-v1.1.csv')
    s = search_from_database(s, database, ppm = 500)