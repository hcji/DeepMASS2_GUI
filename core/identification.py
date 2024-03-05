# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 10:34:38 2023

@author: DELL
"""


import base64

import matchms.filtering as msfilters
import numpy as np
import pandas as pd
from matchms import calculate_scores
from matchms.similarity import CosineGreedy
from molmass import Formula
from rdkit import Chem
from rdkit.Chem import DataStructs, AllChem
from sklearn.metrics.pairwise import cosine_similarity
from spec2vec import SpectrumDocument
from spec2vec.vector_operations import calc_vector

from core.pubchem import retrieve_by_formula, retrieve_by_exact_mass
from core.pubchem import retrieve_by_formula_database, retrieve_by_exact_mass_database
from core.pycdk import IsotopeFromString, IsotopeSimilarity


def spectrum_processing(s):
    """This is how one would typically design a desired pre- and post-
    processing pipeline."""
    s = msfilters.default_filters(s)
    if ("adduct_type" in s.metadata.keys()) and ("adduct" not in s.metadata.keys()):
        s.set("adduct", s.get("adduct_type"))
    s = msfilters.correct_charge(s)
    s = msfilters.add_parent_mass(s)
    s = msfilters.add_losses(s)
    s = msfilters.normalize_intensities(s)
    s = msfilters.select_by_mz(s, mz_from=0, mz_to=1000)
    return s


def get_formula_mass(formula):
    formula = formula.replace("+", "").replace("-", "")
    f = Formula(formula)
    return f.isotope.mass


def search_candidates(s, database, ms1_tolerence=20):
    if "formula" in s.metadata.keys():
        formula = s.metadata["formula"]
        candidate = retrieve_by_formula_database(formula, database)
        return candidate
        if len(candidate) == 0:
            try:
                candidate = retrieve_by_formula(formula)
                return candidate
            except:
                return None

    elif "parent_mass" in s.metadata.keys():
        mass = float(s.metadata["parent_mass"])
        candidate = retrieve_by_exact_mass_database(mass, database, ppm=ms1_tolerence)
        return candidate
        if len(candidate) == 0:
            try:
                candidate = retrieve_by_exact_mass(mass)
                return candidate
            except:
                return None
    else:
        return None


def calc_isotope_score(s, candidate_mol):
    isotope_score = []
    if s.get("isotope_mz") and s.get("isotope_intensity"):
        isotope_mz = (
            base64.b64decode(s.get("isotope_mz").split("'")[1])
            .decode("ascii")
            .replace("\n", "")
        )
        isotope_intensity = (
            base64.b64decode(s.get("isotope_intensity").split("'")[1])
            .decode("ascii")
            .replace("\n", "")
        )
        isotope_mz = [
            float(s)
            for s in isotope_mz.replace("[", "").replace("]", "").split(" ")
            if s != ""
        ]
        isotope_intensity = [
            float(s)
            for s in isotope_intensity.replace("[", "").replace("]", "").split(" ")
            if s != ""
        ]
        isotope_mz = np.array(isotope_mz)

        adduct_mz = isotope_mz[np.argmax(isotope_intensity)] - s.metadata["parent_mass"]
        isotope_mz = isotope_mz - adduct_mz

        isotope_intensity = np.array(isotope_intensity)
        isotope_intensity = isotope_intensity / max(isotope_intensity)
        isotope_pattern = np.vstack((isotope_mz, isotope_intensity)).T
    else:
        return None
    for i in range(len(candidate_mol)):
        try:
            formula = AllChem.CalcMolFormula(candidate_mol[i])
            isotope_ref = IsotopeFromString(formula, minI=0.001)
        except:
            isotope_score.append(0)
        isotope_score.append(IsotopeSimilarity(isotope_pattern, isotope_ref, 10))
    return isotope_score


def calc_deepmass_score(
    s, candidate_mol, reference_mol, query_vector, reference_vector
):
    get_fp = lambda x: AllChem.GetMorganFingerprintAsBitVect(x, radius=2)
    get_sim = lambda x, y: DataStructs.FingerprintSimilarity(x, y)
    get_corr = lambda x, y: cosine_similarity([x], [y])[0][0]

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
    deepmass_score /= np.max(deepmass_score)
    return deepmass_score


def calc_wt_score(s, candidate_mol):
    candidate_mass = [AllChem.CalcExactMolWt(m) for m in candidate_mol]
    diff_mass = np.array([abs(m - float(s.get("parent_mass"))) for m in candidate_mass])
    wt_score = 1 - 50000 * diff_mass / float(s.get("parent_mass"))
    return wt_score


def identify_unknown(s, p, model, references, database):
    candidate = search_candidates(s, database)
    if candidate is None:
        return s
    if len(candidate) == 0:
        return s
    candidate_mol = [Chem.MolFromSmiles(s) for s in candidate["CanonicalSMILES"]]
    query_vector = calc_vector(model, SpectrumDocument(s, n_decimals=2))

    xq = np.array(query_vector).astype("float32")
    I, D = p.knn_query(xq, 300)

    reference_spectrum = np.array(references)[I[0, :]]
    reference_smile = [s.metadata["smiles"] for s in reference_spectrum]
    reference_mol = [Chem.MolFromSmiles(s) for s in reference_smile]
    reference_vector = np.array(p.get_items(I[0, :]))

    candidate_deepmass_score = calc_deepmass_score(
        s, candidate_mol, reference_mol, query_vector, reference_vector
    )
    candidate["DeepMass Score"] = np.round(candidate_deepmass_score, 4)

    if s.get("formula") is None:
        candidate_wt_score = calc_wt_score(s, candidate_mol)
        candidate["MolWt Score"] = np.round(candidate_wt_score, 4)
        if s.get("isotope_mz") and s.get("isotope_intensity"):
            candidate_isotopic_score = calc_isotope_score(s, candidate_mol)
            candidate["Isotope Score"] = np.round(candidate_isotopic_score, 4)
            candidate["Consensus Score"] = (
                0.8 * candidate["DeepMass Score"]
                + 0.1 * candidate["Isotope Score"]
                + 0.1 * candidate["MolWt Score"]
            )
        else:
            candidate["Consensus Score"] = (
                0.8 * candidate["DeepMass Score"] + 0.2 * candidate["MolWt Score"]
            )
        candidate = candidate.sort_values(
            "Consensus Score", ignore_index=True, ascending=False
        )
    else:
        candidate = candidate.sort_values(
            "DeepMass Score", ignore_index=True, ascending=False
        )

    s.set("annotation", candidate)
    s.set("reference", reference_spectrum)
    return s


def match_spectrum(s, precursors, references):
    precursor = s.get("precursor_mz")
    if precursor is None:
        return s
    lb, ub = precursor - 0.05, precursor + 0.05
    li = np.searchsorted(precursors, lb)
    ui = np.searchsorted(precursors, ub)
    if ui <= li:
        return s
    match_scores = calculate_scores(
        references=references[li:ui], queries=[s], similarity_function=CosineGreedy()
    )
    # print(match_scores.scores)
    match_scores = np.array([s[0].tolist()[0] for s in match_scores.scores])
    w = np.argsort(-match_scores)
    match_scores = match_scores[w]
    reference = np.array(references)[li:ui][w]

    annotation, inchikeys = [], []
    for i, r in enumerate(reference):
        mol = Chem.MolFromSmiles(r.get("smiles"))
        score = match_scores[i]
        if mol is None:
            continue
        inchikey = r.get("inchikey")
        if inchikey == "":
            continue
        title = r.get("compound_name")
        smiles = Chem.MolToSmiles(mol)
        try:
            formula = AllChem.CalcMolFormula(mol)
        except:
            formula = ""
        if inchikey not in inchikeys:
            inchikeys.append(inchikey)
            annotation.append([title, formula, smiles, inchikey, score])
    annotation = pd.DataFrame(
        annotation,
        columns=[
            "Title",
            "MolecularFormula",
            "CanonicalSMILES",
            "InChIKey",
            "Matching Score",
        ],
    )

    if s.get("formula") is not None:
        annotation = annotation[annotation["MolecularFormula"] == s.get("formula")]
        annotation = annotation.reset_index(drop=True)
    s.set("annotation", annotation)
    s.set("reference", reference)
    return s


if __name__ == "__main__":
    """
    import hnswlib
    import pickle
    import pandas as pd
    from matchms.importing import load_from_mgf
    from gensim.models import Word2Vec

    model = Word2Vec.load("model/Ms2Vec_allGNPSpositive.hdf5")
    p = hnswlib.Index(space='l2', dim=300)
    p.load_index('data/references_index_positive_spec2vec.bin')
    with open('data/references_spectrums_positive.pickle', 'rb') as file:
        references = pickle.load(file)
    references = np.array(references)
    precursors = [s.get('precursor_mz') for s in references]
    precursors = np.array(precursors)

    spectrums = [s for s in load_from_mgf("D:/DeepMASS2_Data_Processing/Example/CASMI/all_casmi.mgf")]
    s = spectrums[200]
    s = identify_unknown(s, p, model, reference, database)
    """
