import logging
import os
import pickle

import numpy as np
import pandas as pd
from gensim.models import Word2Vec
from hnswlib import Index

from backend.utils.theme import Seafoam
from core.identification import identify_unknown, match_spectrum

# matplotlib.use('Agg')

os.environ["MPLCONFIGDIR"] = os.getcwd() + "/configs/"


seafoam = Seafoam()

default_index_positive = "data/references_index_positive_spec2vec.bin"
default_index_negative = "data/references_index_negative_spec2vec.bin"
default_reference_positive = "data/references_spectrums_positive.pickle"
default_reference_negative = "data/references_spectrums_negative.pickle"
print("Start Loading database")
default_database = pd.read_csv("data/database.csv")
print("Start Loading Word2Vec")
deepmass_positive = Word2Vec.load("model/Ms2Vec_allGNPSpositive.hdf5")
deepmass_negative = Word2Vec.load("model/Ms2Vec_allGNPSnegative.hdf5")
print("Start Loading negative reference")

with open(default_reference_negative, "rb") as file:
    reference_negative = pickle.load(file)
print("Start Loading positive reference")
with open(default_reference_positive, "rb") as file:
    reference_positive = pickle.load(file)
print("Start Loading hnsw index")
index_negative = Index(space="l2", dim=300)
index_negative.load_index(default_index_negative)

index_positive = Index(space="l2", dim=300)
index_positive.load_index(default_index_positive)

precursors_positive = np.array([s.get("precursor_mz") for s in reference_positive])
precursors_negative = np.array([s.get("precursor_mz") for s in reference_negative])
print("Finish!!!")
print("-" * 100)


def identify_pos(spectrum):
    return identify_unknown(
        spectrum,
        index_positive,
        deepmass_positive,
        reference_positive,
        default_database,
    )


def identify_neg(spectrum):
    return identify_unknown(
        spectrum,
        index_negative,
        deepmass_negative,
        reference_negative,
        default_database,
    )


def match_pos(spectrum):
    return match_spectrum(spectrum, precursors_positive, reference_positive)


def match_neg(spectrum):
    return match_spectrum(spectrum, precursors_negative, reference_negative)


def id_spectrum_list(spectrum_list, progress=None, is_deepmass=True):
    res = []
    if is_deepmass:
        for s in progress.tqdm(spectrum_list):
            logging.info(f"")
            sn = None
            if "ionmode" in s.metadata.keys():
                if s.metadata["ionmode"] == "negative":
                    sn = identify_neg(s)
                else:
                    sn = identify_pos(s)
            else:
                sn = identify_pos(s)
            res.append(sn)
    else:
        for s in progress.tqdm(spectrum_list):
            sn = None
            if "ionmode" in s.metadata.keys():
                if s.metadata["ionmode"] == "negative":
                    sn = match_neg(s)
                else:
                    sn = match_pos(s)
            else:
                sn = match_pos(s)
            res.append(sn)
    return res
