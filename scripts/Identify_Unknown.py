# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 08:41:08 2024

@author: DELL
"""


import numpy as np
import pandas as pd
from matchms.Spectrum import Spectrum
from matchms.filtering import default_filters, normalize_intensities

data_GH11 = np.array(pd.read_csv("D:/DeepMASS2_Data_Processing/Example/Unknown/GH11.csv", header=None))
spec_GH11 = Spectrum(mz = data_GH11[:,0], intensities = data_GH11[:,1],
                  metadata = {'name':'GH11', 'precursor_mz': 287.1751, 'formula': 'C17H22N2O2', 'parent_mass': 0})
spec_GH11 = normalize_intensities(default_filters(spec_GH11))

data_G7H = np.array(pd.read_csv("D:/DeepMASS2_Data_Processing/Example/Unknown/G7H.csv", header=None))
spec_G7H = Spectrum(mz = data_G7H[:,0], intensities = data_G7H[:,1],
                  metadata = {'name':'G7H', 'precursor_mz': 293.1646, 'formula': 'C19H12N2O', 'parent_mass': 0})
spec_G7H = normalize_intensities(default_filters(spec_G7H))

data_GH17 = np.array(pd.read_csv("D:/DeepMASS2_Data_Processing/Example/Unknown/GH17.csv", header=None))
spec_GH17 = Spectrum(mz = data_GH17[:,0], intensities = data_GH17[:,1],
                  metadata = {'name':'GH17', 'precursor_mz': 231.1490, 'formula': 'C14H18N2O', 'parent_mass': 0})
spec_GH17 = normalize_intensities(default_filters(spec_GH17))


import hnswlib
import pickle
from gensim.models import Word2Vec
from rdkit import Chem
from rdkit.Chem import Draw
from core.annotating.candidates import search_from_database, search_from_pubchem
from core.annotating.structure import calc_deepmass_score

database = pd.read_csv('data/DeepMassStructureDB-v1.1.csv')
model = Word2Vec.load("model/Ms2Vec_allGNPSpositive.hdf5")
p = hnswlib.Index(space='l2', dim=300) 
p.load_index('data/references_index_positive_spec2vec.bin')
with open('data/references_spectrums_positive.pickle', 'rb') as file:
    references = pickle.load(file)


# G7H
spec_G7H = search_from_pubchem(spec_G7H)
deepmass_score, reference_spectrum = calc_deepmass_score(spec_G7H, p, model, references)

true_annotation = pd.Series(['G7H', 'C19H12N2O', 'C1=CC=CC2=C1C=C([N]2CC3=CC=CC=C3)N4CCOCC4', 'YARLJYDZSKQYIM-UHFFFAOYSA-N', np.nan, np.nan, np.nan], 
                             index = spec_G7H.metadata['annotation'].columns)
spec_G7H.set('annotation', spec_G7H.metadata['annotation'].append(true_annotation, ignore_index = True))
deepmass_score_new, reference_spectrum = calc_deepmass_score(spec_G7H, p, model, references)

Draw.MolToImage(Chem.MolFromSmiles('C1=CC=CC2=C1C=C([N]2CC3=CC=CC=C3)N4CCOCC4'))
Draw.MolToImage(Chem.MolFromSmiles(reference_spectrum[0].get('smiles')))

# GH11
spec_GH11 = search_from_database(spec_GH11, database)
deepmass_score, reference_spectrum = calc_deepmass_score(spec_GH11, p, model, references)

true_annotation = pd.Series(['GH11', 'C17H22N2O2', 'C1=CC=CC3=C1C=C(N2CCC(CC2)C(=O)OCC)[N]3C', 'HFDRVGLLZYKILL-UHFFFAOYSA-N', np.nan, np.nan, np.nan, np.nan], 
                             index = spec_GH11.metadata['annotation'].columns)
spec_GH11.set('annotation', spec_GH11.metadata['annotation'].append(true_annotation, ignore_index = True))
deepmass_score_new, reference_spectrum = calc_deepmass_score(spec_GH11, p, model, references)

Draw.MolToImage(Chem.MolFromSmiles('C1=CC=CC3=C1C=C(N2CCC(CC2)C(=O)OCC)[N]3C'))
Draw.MolToImage(Chem.MolFromSmiles(reference_spectrum[0].get('smiles')))


# GH17
spec_GH17 = search_from_database(spec_GH17, database)
deepmass_score, reference_spectrum = calc_deepmass_score(spec_GH17, p, model, references)

true_annotation = pd.Series(['GH17', 'C14H18N2O', 'C1(=CC=CC3=C1C=C(N2CCOCC2)[N]3C)C', 'UEFBWULSZPHZBE-UHFFFAOYSA-N', np.nan, np.nan, np.nan, np.nan], 
                             index = spec_GH17.metadata['annotation'].columns)
spec_GH17.set('annotation', spec_GH17.metadata['annotation'].append(true_annotation, ignore_index = True))
deepmass_score_new, reference_spectrum = calc_deepmass_score(spec_GH17, p, model, references)

Draw.MolToImage(Chem.MolFromSmiles('C1(=CC=CC3=C1C=C(N2CCOCC2)[N]3C)C'))
Draw.MolToImage(Chem.MolFromSmiles(reference_spectrum[5].get('smiles')))

