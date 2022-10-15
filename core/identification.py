# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 11:27:02 2022

@author: DELL
"""


import numpy as np
from rdkit import Chem
from rdkit.Chem import DataStructs, AllChem
from sklearn.metrics.pairwise import cosine_similarity
import matchms.filtering as msfilters

from core.pubchem import retrieve_by_formula, retrieve_by_exact_mass
from core.pubchem import retrieve_by_formula_database, retrieve_by_exact_mass_database


def spectrum_processing(s):
    """This is how one would typically design a desired pre- and post-
    processing pipeline."""
    s = msfilters.default_filters(s)
    s = msfilters.add_parent_mass(s)
    s = msfilters.normalize_intensities(s)
    s = msfilters.select_by_mz(s, mz_from=0, mz_to=1000)
    return s


def identify_unknown(s, p, n, database, priority, model, reference, chemical_space, in_silicon_only=False):
    """
    Example:
        import hnswlib
        import pandas as pd
        from ms2deepscore import MS2DeepScore
        from ms2deepscore.models import load_model
        from matchms.importing import load_from_mgf
        
        n = 30
        priority = ['HMDB', 'KNApSAcK', 'BLEXP']
        spectrums = [s for s in load_from_mgf('example/casmi_2022_priority_new.mgf')]
        model = load_model("model/MS2DeepScore_allGNPSpositive.hdf5")
        model = MS2DeepScore(model)
        database = pd.read_csv('data/MsfinderStructureDB-VS15-plus-GNPS.csv')
        p = hnswlib.Index(space='l2', dim=200) 
        p.load_index('data/references_index_positive.bin')
        reference = np.load('data/references_spectrums_positive.npy', allow_pickle=True)
        s = spectrums[2]
        sn = identify_unknown(s, p, n, database, priority, model, reference, 'custom', True)
    """
    s = spectrum_processing(s)
    query_vector = model.calculate_vectors([s])
    xq = np.array(query_vector).astype('float32')
    I, D = p.knn_query(xq, n)
    
    get_fp = lambda x: AllChem.GetMorganFingerprintAsBitVect(x, radius=2)
    get_sim = lambda x, y: DataStructs.FingerprintSimilarity(x, y)
    get_corr = lambda x, y: cosine_similarity([x], [y])[0][0]
    
    if 'formula' in s.metadata.keys():
        formula = s.metadata['formula']
        if chemical_space == 'biodatabase':
            candidate = retrieve_by_formula_database(formula, database, priority = priority)
        elif chemical_space == 'biodatabase plus':
            candidate = retrieve_by_formula_database(formula, database, priority = priority)
            if len(candidate) == 0:
                try:
                    candidate = retrieve_by_formula(formula)
                except:
                    candidate = retrieve_by_formula_database(formula, database, priority = [])
        elif chemical_space == 'custom':
            candidate = retrieve_by_formula_database(formula, database, priority = [])
        else:
            try:
                candidate = retrieve_by_formula(formula)
            except:
                candidate = retrieve_by_formula_database(formula, database, priority = [])
             
    elif 'parent_mass' in s.metadata.keys():
        mass = s.metadata['parent_mass']
        if chemical_space == 'biodatabase':
            candidate = retrieve_by_exact_mass_database(mass, database, ppm = 10, priority = priority)
        elif chemical_space == 'biodatabase plus':
            candidate = retrieve_by_exact_mass_database(mass, database, ppm = 10, priority = priority)
            if len(candidate) == 0:
                try:
                    candidate = retrieve_by_exact_mass(mass)
                except:
                    candidate = retrieve_by_exact_mass_database(mass, database, ppm = 10, priority = [])        
        elif chemical_space == 'custom':
            candidate = retrieve_by_exact_mass_database(mass, database, ppm = 10, priority = [])
        else:
            try:
                candidate = retrieve_by_exact_mass(mass)
            except:
                candidate = retrieve_by_exact_mass_database(mass, database, ppm = 10, priority = [])
    else:
        return s
    
    if len(candidate) == 0:
        return s
    reference_spectrum = np.array(reference)[I[0,:]]
    reference_smile = [s.metadata['smiles'] for s in reference_spectrum]
    reference_mol = [Chem.MolFromSmiles(s) for s in reference_smile]
    k, reference_fp = [], []
    for i, m in enumerate(reference_mol):
        if in_silicon_only:
            if 'inchikey' in s.metadata.keys():
                que_key = s.metadata['inchikey']
                ref_key = reference_spectrum[i].metadata['inchikey'][:14]
                if que_key == ref_key:
                    continue
        try:
            reference_fp.append(get_fp(m))
            k.append(i)
        except:
            pass
    reference_spectrum = np.array(reference_spectrum)[np.array(k)]
    reference_smile = np.array(reference_smile)[np.array(k)]
    reference_vec = p.get_items(I[0,:])
    reference_corr = [get_corr(query_vector[0,:], v) for v in reference_vec]
    
    candidate_mol = [Chem.MolFromSmiles(s) for s in candidate['CanonicalSMILES']]
    k, candidate_fp = [], []
    for i, m in enumerate(candidate_mol):
        try:
            candidate_fp.append(get_fp(m))
            k.append(i)
        except:
            pass
    if len(k) == 0:
        return s
    candidate = candidate.loc[np.array(k),:].reset_index(drop = True)
    candidate_fp_sim = np.array([[get_sim(f1, f2) for f1 in reference_fp] for f2 in candidate_fp])
    candidate_fp_score = [np.sqrt(reference_corr * s) for s in candidate_fp_sim]
    candidate_fp_deepmass = np.sum(candidate_fp_score, axis = 1)
    candidate['DeepMass Score'] = candidate_fp_deepmass / n
    candidate = candidate.sort_values('DeepMass Score', ignore_index = True, ascending = False)
    candidate['DeepMass Score'] = np.round(candidate['DeepMass Score'], 4)
    
    s.set('reference', reference_spectrum)
    s.set('annotation', candidate)
    return s
    
    

