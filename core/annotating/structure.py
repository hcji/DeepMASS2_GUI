# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 08:32:04 2023

@author: DELL
"""


import numpy as np
import pandas as pd
from matchms import calculate_scores
from matchms.similarity import CosineGreedy
from rdkit import Chem
from rdkit.Chem import DataStructs, AllChem
from sklearn.metrics.pairwise import cosine_similarity
from spec2vec import SpectrumDocument
from spec2vec.vector_operations import calc_vector


def check_inputs(s):
    if s.get("annotation") is None:
        return False
    elif s.get("parent_mass") is None:
        return False
    elif s.get("precursor_mz") is None:
        return False
    else:
        return True


def calc_deepmass_score(s, p, model, references):
    if not check_inputs(s):
        return None
    else:
        pass
    get_fp = lambda x: AllChem.GetMorganFingerprintAsBitVect(x, radius=2)
    get_sim = lambda x, y: DataStructs.FingerprintSimilarity(x, y)
    get_corr = lambda x, y: cosine_similarity([x], [y])[0][0]
    candidate_mol = [
        Chem.MolFromSmiles(smi) for smi in s.get("annotation")["CanonicalSMILES"]
    ]

    query_vector = calc_vector(
        model, SpectrumDocument(s, n_decimals=2), allowed_missing_percentage=100
    )
    xq = np.array(query_vector).astype("float32")
    I, D = p.knn_query(xq, 300)
    reference_spectrum = np.array(references)[I[0, :]]
    reference_spectrum = [s for s in reference_spectrum if s.get("smiles") is not None]
    reference_smile = [s.metadata["smiles"] for s in reference_spectrum]
    reference_mol = [Chem.MolFromSmiles(s) for s in reference_smile]
    reference_vector = np.array(p.get_items(I[0, :]))

    k, reference_fp = [], []
    for i, m in enumerate(reference_mol):
        try:
            reference_fp.append(get_fp(m))
            k.append(i)
        except:
            pass
    if len(k) != len(reference_mol):
        k = np.array(k)
        reference_mol = np.array(reference_mol)[k]
        reference_vector = np.array(reference_vector)[k, :]

    if len(candidate_mol) == 0:
        return None

    deepmass_score = []
    for i in range(len(candidate_mol)):
        try:
            candidate_fp_i = get_fp(candidate_mol[i])
        except:
            deepmass_score.append(0)
        candidate_vecsim_i = [
            get_corr(query_vector, reference_vector_)
            for reference_vector_ in reference_vector
        ]
        candidate_vecsim_i = np.array(candidate_vecsim_i)
        candidate_fpsim_i = [
            get_sim(candidate_fp_i, reference_fp_) for reference_fp_ in reference_fp
        ]
        candidate_fpsim_i = np.array(candidate_fpsim_i)
        top20 = np.argsort(-np.array(candidate_fpsim_i))[:20]
        candidate_score_i = np.sqrt(
            np.sum(candidate_vecsim_i[top20] * candidate_fpsim_i[top20])
        )
        deepmass_score.append(candidate_score_i / 20)
    deepmass_score = np.array(deepmass_score)
    deepmass_score /= np.max(deepmass_score) + 10**-10
    deepmass_score = dict(zip(s.get("annotation")["InChIKey"], deepmass_score))
    return deepmass_score, reference_spectrum


def calc_matchms_score(s, precursors, references):
    if not check_inputs(s):
        return None
    else:
        pass
    precursor = s.get("precursor_mz")
    lb, ub = precursor - 0.05, precursor + 0.05
    li = np.searchsorted(precursors, lb)
    ui = np.searchsorted(precursors, ub)

    reference_spectrum = np.array(references)[li:ui]
    reference_spectrum = np.array(
        [s for s in reference_spectrum if s.get("inchikey") is not None]
    )
    if len(reference_spectrum) == 0:
        return None, None

    match_scores = calculate_scores(
        references=reference_spectrum, queries=[s], similarity_function=CosineGreedy()
    )
    # match_scores = match_scores.scores.to_array()
    match_scores = np.array([s[0].tolist()[0] for s in match_scores.scores.to_array()])
    w = np.argsort(-match_scores)
    match_scores = match_scores[w]
    reference_spectrum = reference_spectrum[w]
    reference_inchikey = np.array([s.get("inchikey")[:14] for s in reference_spectrum])

    candidate_inchiky = s.get("annotation")["InChIKey"]
    candidate_scores = []
    for k in candidate_inchiky:
        k = str(k)[:14]
        w = np.where(reference_inchikey == k)[0]
        if len(w) == 0:
            candidate_scores.append(-float("inf"))
        else:
            candidate_scores.append(np.max(match_scores[w]))
    candidate_scores = dict(zip(candidate_inchiky, candidate_scores))
    return candidate_scores, reference_spectrum


def calc_multiomics_score(s, association_data):
    associated_id = s.get("associated_gene")
    if associated_id is None:
        return None
    else:
        get_fp = lambda x: AllChem.GetMorganFingerprintAsBitVect(x, radius=2)
        get_sim = lambda x, y: DataStructs.FingerprintSimilarity(x, y)
        associated_id = associated_id.split(",")

    associated_smiles = []
    for gene in associated_id:
        j = np.where(association_data["Gene Name"] == gene)[0]
        if len(j) > 0:
            for jj in j:
                associated_smiles += list(
                    association_data.loc[jj, "Association SMILES"]
                )
    associated_smiles = np.array(associated_smiles)

    if s.get("smiles") is not None:
        try:
            true_molwt = AllChem.CalcExactMolWt(Chem.MolFromSmiles(s.get("smiles")))
            associated_molwt = np.array(
                [
                    AllChem.CalcExactMolWt(Chem.MolFromSmiles(s))
                    for s in associated_smiles
                ]
            )
            k = np.where(associated_molwt != true_molwt)[0]
            associated_smiles = associated_smiles[k]
        except:
            pass
    if len(associated_smiles) == 0:
        return None

    k, associated_fp = [], []
    for i, smi in enumerate(associated_smiles):
        try:
            mol = Chem.MolFromSmiles(smi)
            associated_fp.append(get_fp(mol))
            k.append(i)
        except:
            pass

    multiomics_score = []
    candidate_mol = [
        Chem.MolFromSmiles(smi) for smi in s.get("annotation")["CanonicalSMILES"]
    ]
    for i in range(len(candidate_mol)):
        try:
            candidate_fp_i = get_fp(candidate_mol[i])
        except:
            multiomics_score.append(0)

        candidate_fpsim_i = [
            get_sim(candidate_fp_i, associated_fp_) for associated_fp_ in associated_fp
        ]
        candidate_fpsim_i = np.array(candidate_fpsim_i)

        top5 = np.argsort(-np.array(candidate_fpsim_i))[:5]
        candidate_score_i = np.sum(candidate_fpsim_i[top5])
        multiomics_score.append(candidate_score_i / 5)
    multiomics_score = np.array(multiomics_score)
    multiomics_score /= np.max(multiomics_score)
    multiomics_score = dict(zip(s.get("annotation")["InChIKey"], multiomics_score))
    return multiomics_score


if __name__ == "__main__":
    import hnswlib
    import pickle
    from gensim.models import Word2Vec
    from core.importing.load_from_files import load_from_files
    from core.annotating.candidates import search_from_database

    # s = load_from_files(['example/minimum_example.mat'])[0]
    s = load_from_files(["example/all_casmi.mgf"])[2]
    database = pd.read_csv("data/DeepMassStructureDB-v1.1.csv")
    s = search_from_database(s, database, ppm=500)

    model = Word2Vec.load("model/Ms2Vec_allGNPSpositive.hdf5")
    p = hnswlib.Index(space="l2", dim=300)
    p.load_index("data/references_index_positive_spec2vec.bin")
    with open("data/references_spectrums_positive.pickle", "rb") as file:
        references = pickle.load(file)
    references = np.array(references)

    precursors = [s.get("precursor_mz") for s in references]
    precursors = np.array(precursors)

    association_data = pd.read_pickle("data/metabolite_gene_associations.pickle")

    calc_deepmass_score(s, p, model, references)
    calc_matchms_score(s, precursors, references)
