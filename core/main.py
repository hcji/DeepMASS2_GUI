# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 14:24:31 2023

@author: DELL
"""


from core.annotating.candidates import search_from_database, search_from_pubchem
from core.annotating.formula import calc_formula_score
from core.annotating.structure import calc_deepmass_score, calc_matchms_score


def identify_unknown(s, p, model, references, database):
    s = search_from_database(s, database, ppm = 50)
    candidate = s.get('annotation')
    if candidate is None:
        try:
            s = search_from_pubchem(s, ppm = 50)
        except:
            return s
    
    candidate = s.get('annotation')
    if candidate is None:
        return s
    
    formula_score = calc_formula_score(s)
    structure_score, reference_spectrum = calc_deepmass_score(s, p, model, references)
    for i in candidate.index:
        k = candidate.loc[i, 'InChIKey']
        f = candidate.loc[i, 'MolecularFormula']
        candidate.loc[i, 'Formula Score'] = formula_score[f]
        candidate.loc[i, 'Structure Score'] = structure_score[k]
        candidate.loc[i, 'Consensus Score'] = 0.3*formula_score[f] + 0.7*structure_score[k]
    candidate = candidate.sort_values('Consensus Score', ignore_index = True, ascending = False)
    s.set('annotation', candidate)
    s.set('reference', reference_spectrum)
    return s


def match_spectrum(s, precursors, references, database):
    s = search_from_database(s, database, ppm = 50)
    candidate = s.get('annotation')
    if candidate is None:
        s = search_from_pubchem(s, ppm = 500)
    
    candidate = s.get('annotation')
    if candidate is None:
        return s
    if len(candidate) == 0:
        return s
    
    formula_score = calc_formula_score(s)
    structure_score, reference_spectrum = calc_matchms_score(s, precursors, references)
    
    if structure_score is None:
        return s
    
    for i in candidate.index:
        k = candidate.loc[i, 'InChIKey']
        f = candidate.loc[i, 'MolecularFormula']
        candidate.loc[i, 'Formula Score'] = formula_score[f]
        candidate.loc[i, 'Structure Score'] = structure_score[k]
        candidate.loc[i, 'Consensus Score'] = 0.3*formula_score[f] + 0.7*structure_score[k]
    candidate = candidate.sort_values('Consensus Score', ignore_index = True, ascending = False)
    s.set('annotation', candidate)
    s.set('reference', reference_spectrum)    
    return s


if __name__ == '__main__':

    '''
    import hnswlib
    import pickle
    import numpy as np
    import pandas as pd
    from gensim.models import Word2Vec
    from core.importing.load_from_files import load_from_files
    
    s = load_from_files(['example/minimum_example.mat'])[0]
    database = pd.read_csv('data/DeepMassStructureDB-v1.1.csv')
    s = search_from_database(s, database, ppm = 500)
    model = Word2Vec.load("model/Ms2Vec_allGNPSpositive.hdf5")
    p = hnswlib.Index(space='l2', dim=300) 
    p.load_index('data/references_index_positive_spec2vec.bin')
    with open('data/references_spectrums_positive.pickle', 'rb') as file:
        references = pickle.load(file)
    references = np.array(references)

    precursors = [s.get('precursor_mz') for s in references]
    precursors = np.array(precursors)
    
    s1 = identify_unknown(s, p, model, references, database)
    s2 = match_spectrum(s, precursors, references, database)
    '''