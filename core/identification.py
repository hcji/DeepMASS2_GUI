# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 11:27:02 2022

@author: DELL
"""


import base64
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs, AllChem
from sklearn.metrics.pairwise import cosine_similarity

import matchms.filtering as msfilters
from molmass import Formula
from spec2vec import SpectrumDocument
from spec2vec.vector_operations import calc_vector
from ms2deepscore.MS2DeepScore import MS2DeepScore
from gensim.models.word2vec import Word2Vec

from pycdk.pycdk import IsotopeFromString, IsotopeSimilarity
from core.pubchem import retrieve_by_formula, retrieve_by_exact_mass
from core.pubchem import retrieve_by_formula_database, retrieve_by_exact_mass_database

structure_db = pd.read_csv('data/DeepMassStructureDB-v1.0.csv')

def spectrum_processing(s):
    """This is how one would typically design a desired pre- and post-
    processing pipeline."""
    s = msfilters.default_filters(s)
    if ('adduct_type' in s.metadata.keys()) and ('adduct' not in s.metadata.keys()):
        s.set('adduct', s.get('adduct_type'))
    s = msfilters.correct_charge(s)
    s = msfilters.add_parent_mass(s)
    s = msfilters.add_losses(s)
    s = msfilters.normalize_intensities(s)
    s = msfilters.select_by_mz(s, mz_from=0, mz_to=1000)
    return s


def get_formula_mass(formula):
    f = Formula(formula)
    return f.isotope.mass


def calc_isotope_score(s, formula):
    # print(formula)
    if s.get('isotope_mz') and s.get('isotope_intensity'):
        isotope_mz = base64.b64decode(s.get('isotope_mz').split("'")[1]).decode("ascii").replace('\n', '')
        isotope_intensity = base64.b64decode(s.get('isotope_intensity').split("'")[1]).decode("ascii").replace('\n', '')
        isotope_mz = [float(s) for s in isotope_mz.replace('[', '').replace(']', '').split(' ') if s != '']
        isotope_intensity = [float(s) for s in isotope_intensity.replace('[', '').replace(']', '').split(' ') if s != '']
        isotope_mz = np.array(isotope_mz)

        adduct_mz = isotope_mz[np.argmax(isotope_intensity)] - s.metadata['parent_mass']
        isotope_mz = isotope_mz - adduct_mz
        
        isotope_intensity = np.array(isotope_intensity)
        isotope_intensity = isotope_intensity / max(isotope_intensity)
        isotope_pattern = np.vstack((isotope_mz, isotope_intensity)).T
    else:
        return 0
    try:
        isotope_ref = IsotopeFromString(formula, minI=0.001)
    except:
        return 0
    return IsotopeSimilarity(isotope_pattern, isotope_ref, 10)


def identify_unknown(s, p, n_ref, n_neb, database, priority, model, reference, chemical_space, in_silicon_only=True, ms1_tolerence = 10):
    """
    Example:
        import hnswlib
        import pickle
        import pandas as pd
        from ms2deepscore import MS2DeepScore
        from ms2deepscore.models import load_model
        from matchms.importing import load_from_mgf
        
        n_ref = 300
        n_neb = 20
        priority = []
        spectrums = [s for s in load_from_mgf('example/Plasma/ms_ms_plasma.mgf')]
        
        model = Word2Vec.load("model/Ms2Vec_allGNPSnegative.hdf5_iter_30.model")
        database = pd.read_csv('data/MsfinderStructureDB-VS15-plus-GNPS.csv')
        p = hnswlib.Index(space='l2', dim=300) 
        p.load_index('data/references_index_negative_spec2vec.bin')
        with open('data/references_spectrums_negative.pickle', 'rb') as file:
            reference = pickle.load(file)
        s = spectrums[427]
        sn = identify_unknown(s, p, n_ref, n_neb, database, priority, model, reference, 'biodatabase', True)
    """
    s = spectrum_processing(s)
    try:
        if type(model) == MS2DeepScore:
            query_vector = model.calculate_vectors([s])[0,:]
        elif type(model) == Word2Vec:
            query_vector = calc_vector(model, SpectrumDocument(s, n_decimals=2))
        else:
            return s
        s.set('query_vector', list(query_vector))
    except:
        s.set('query_vector', list())
        return s
    xq = np.array(query_vector).astype('float32')
    I, D = p.knn_query(xq, n_ref)
    
    get_fp = lambda x: AllChem.GetMorganFingerprintAsBitVect(x, radius=2)
    get_sim = lambda x, y: DataStructs.FingerprintSimilarity(x, y)
    get_corr = lambda x, y: cosine_similarity([x], [y])[0][0]
    
    if 'formula' in s.metadata.keys():
        formula = s.metadata['formula']
        if chemical_space == 'biodatabase':
            candidate = retrieve_by_formula_database(formula, database)
        elif chemical_space == 'biodatabase plus':
            candidate = retrieve_by_formula_database(formula, database)
            if len(candidate) == 0:
                try:
                    candidate = retrieve_by_formula(formula)
                except:
                    pass
        elif chemical_space == 'custom':
            candidate = retrieve_by_formula_database(formula, database)
        else:
            try:
                candidate = retrieve_by_formula(formula)
            except:
                candidate = retrieve_by_formula_database(formula, database)
             
    elif 'parent_mass' in s.metadata.keys():
        try:
            mass = float(s.metadata['parent_mass'])
        except:
            return s
        if chemical_space == 'biodatabase':
            candidate = retrieve_by_exact_mass_database(mass, database, ppm = ms1_tolerence)
        elif chemical_space == 'biodatabase plus':
            candidate = retrieve_by_exact_mass_database(mass, database, ppm = ms1_tolerence)
            if len(candidate) == 0:
                try:
                    candidate = retrieve_by_exact_mass(mass)
                except:
                    pass
        elif chemical_space == 'custom':
            candidate = retrieve_by_exact_mass_database(mass, database, ppm = ms1_tolerence)
        else:
            try:
                candidate = retrieve_by_exact_mass(mass)
            except:
                candidate = retrieve_by_exact_mass_database(mass, database, ppm = ms1_tolerence)
    else:
        return s
    
    if len(candidate) == 0:
        return s
    
    priority_key = []
    if len(priority) > 0:
        k = set()
        for pr in priority:
            k = k | set(np.where(structure_db[pr].values.astype(str) != 'nan')[0])
        priority_key = list(structure_db.loc[list(k), 'InChIkey'])
    
    reference_spectrum = np.array(reference)[I[0,:]]
    reference_smile = [s.metadata['smiles'] for s in reference_spectrum]
    reference_mol = [Chem.MolFromSmiles(s) for s in reference_smile]
    k, reference_fp = [], []
    for i, m in enumerate(reference_mol):
        if reference_smile[i] == '':
            continue
        if reference_spectrum[i].metadata['ionmode'] != s.metadata['ionmode']:
            continue
        # in-silicon only means reference structures cannot include any candidates
        if in_silicon_only:
            if abs(float(s.metadata['parent_mass']) - reference_spectrum[i].metadata['parent_mass']) < 0.01:
                continue
        try:
            reference_fp.append(get_fp(m))
            k.append(i)
        except:
            pass
    if len(k) == 0:
        return s
    
    reference_spectrum = np.array(reference_spectrum)[np.array(k)]
    reference_smile = np.array(reference_smile)[np.array(k)]
    reference_vec = p.get_items(I[0, np.array(k)])
    reference_corr = [get_corr(query_vector, v) for v in reference_vec]
    
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
    candidate_fp_deepmass = np.array([np.sum(-np.sort(-s)[:n_neb]) for s in candidate_fp_score])
    candidate['DeepMass Score'] = candidate_fp_deepmass / n_neb
    
    formula = candidate['MolecularFormula']
    mass = [get_formula_mass(f) for f in formula]
    if 'parent_mass' in s.metadata.keys():
        diff = np.array([abs(m - s.metadata['parent_mass']) for m in mass])
        candidate_wt_score = 1 - 50000 * diff / s.metadata['parent_mass']
    else:
        candidate_wt_score = np.zeros(len(formula))
    candidate['MolWt Score'] = np.round(candidate_wt_score, 4)
    
    candidate_iso_score = [calc_isotope_score(s, f) for f in candidate['MolecularFormula'].values]
    candidate['Isotope Score'] = np.round(candidate_iso_score, 4)
    
    candidate['Database Score'] = np.array([i in priority_key for i in candidate['InChIKey']]).astype(int)
    
    k = np.argsort(-candidate['DeepMass Score'])
    deepmass_score = np.array([-np.sort(-s) for s in candidate_fp_score])[k,:]
    
    candidate['DeepMass Score'] = np.round(candidate['DeepMass Score'], 4)
    candidate['Consensus Score'] = candidate['DeepMass Score'] + 0.1*candidate['Isotope Score'] + 0.1*candidate['MolWt Score'] + 0.1 * candidate['Database Score']
    candidate['Consensus Score'] /= np.max(candidate['Consensus Score'])
    candidate = candidate.sort_values('Consensus Score', ignore_index = True, ascending = False)
    
    reference_shortkey = [s.get('inchikey')[:14] for s in reference_spectrum]
    candidate_in_reference = [str(s[:14] in reference_shortkey) for s in candidate['InChIKey']]
    candidate['In Reference'] = candidate_in_reference
    
    s.set('reference', reference_spectrum)
    s.set('annotation', candidate)
    s.set('deepmass_score', deepmass_score)
    return s
    