# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 10:16:57 2023

@author: DELL
"""


import numpy as np

from core.annotating.subformula.mass_spectrum import MassSpectrum
from core.pycdk import IsotopeFromString, IsotopeSimilarity, getFormulaExactMass


def check_inputs(s):
    if s.get('annotation') is None:
        return False
    elif s.isotopic_pattern is None:
        return False
    elif s.get('parent_mass') is None:
        return False
    elif s.get('precursor_mz') is None:
        return False
    else:
        return True


def calc_isotopic_score(s):
    if not check_inputs(s):
        return None
    else:
        pass
    formula = set(s.get('annotation').loc[:,'MolecularFormula'])
    adduct_mz = s.get('precursor_mz') - s.get('parent_mass')
    isotope_mz = s.isotopic_pattern.mz
    isotope_mz = isotope_mz - adduct_mz
    isotope_intensity = s.isotopic_pattern.intensities
    isotope_intensity = isotope_intensity / max(isotope_intensity)
    isotope_pattern = np.vstack((isotope_mz, isotope_intensity)).T
    isotope_score = {}
    for f in formula:
        try:
            isotope_ref = IsotopeFromString(f, minI=0.001)
            isotope_score[f] = IsotopeSimilarity(isotope_pattern, isotope_ref, 10)
        except:
            isotope_score[f] = 0
    return isotope_score


def calc_exact_mass_score(s):
    parent_mass = s.get('parent_mass')
    formula = set(s.get('annotation').loc[:,'MolecularFormula'])
    exact_mass_score = {}
    for f in formula:
        try:
            formula_mass = getFormulaExactMass(f)
            diff_mass = abs(formula_mass - parent_mass)
            exact_mass_score[f] = 1 - 1000 * diff_mass / parent_mass
        except:
            exact_mass_score[f] = 0
    return exact_mass_score


def calc_fragmentation_tree_score(s):
    if not check_inputs(s):
        return s
    else:
        pass
    formula = list(set(s.get('annotation').loc[:,'MolecularFormula']))
    formula_mass = np.array([getFormulaExactMass(f) for f in formula])
    formula_mass_ppm = np.max(np.abs(formula_mass - s.get('parent_mass')) / s.get('parent_mass')  * 10**6)
    fragmentation = MassSpectrum(dict(zip(s.mz, s.intensities)), s.get('precursor_mz'), formula_mass_ppm)
    annotations = fragmentation.compute_with_candidate_formula(formula, fragmentation.product_scoring_function)
    scores = [fragmentation.product_scoring_function(s) for s in annotations]
    fragmentation_tree_score = dict(zip(formula, scores))
    return fragmentation_tree_score


def calc_formula_score(s):
    isotope_score = calc_isotopic_score(s)
    exact_mass_score = calc_exact_mass_score(s)
    # fragmentation_tree_score = calc_fragmentation_tree_score
    if isotope_score is None:
        return exact_mass_score
    
    keys = exact_mass_score.keys()
    consesus_score = {}
    for k in keys:
        consesus_score[k] = (exact_mass_score[k] + isotope_score[k]) / 2
    return consesus_score


if __name__ == '__main__':
    
    '''
    from core.importing.load_from_files import load_from_files
    from core.annotating.candidates import search_from_database
    
    s = load_from_files(['example/minimum_example.mat'])[0]
    database = pd.read_csv('data/DeepMassStructureDB-v1.1.csv')
    s = search_from_database(s, database, ppm = 500)
    '''
    