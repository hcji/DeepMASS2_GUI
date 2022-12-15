# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 16:19:45 2022

@author: DELL
"""

import os
import pickle
import numpy as np
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import rdmolfiles, inchi
from matchms.importing import load_from_msp

path_data = os.path.join('D:/All_MSDatabase/NIST2020')

file_mol = os.path.join(path_data, 'hr_msms_nist.MOL')
file_spec = os.path.join(path_data, 'hr_msms_nist.MSP')

spectrums = [s for s in tqdm(load_from_msp(file_spec))]
np.save(os.path.join(path_data, 'preprocessed_spectrums.npy'), spectrums)

def add_mol_info(s):
    if s.get('smiles'):
        return s
    i = s.get('id')
    f = file_mol + '/ID{}.mol'.format(i)
    try:
        m = rdmolfiles.MolFromMolFile(f)
    except:
        return None
    if m is None:
        return None
    smi = Chem.MolToSmiles(m)
    inchikey = inchi.MolToInchiKey(m)
    s = s.set('smiles', smi)
    s = s.set('inchikey', inchikey)
    return s


spectrums = [add_mol_info(s) for s in tqdm(spectrums) if s is not None]
np.save(os.path.join(path_data, 'preprocessed_spectrums.npy'), spectrums)


from matchms.filtering import default_filters
from matchms.filtering import add_parent_mass
from matchms.filtering.load_adducts import load_adducts_dict

adducts_dict = load_adducts_dict()
adducts_dict = {k: v for k, v in adducts_dict.items() if k in ['[M+H]+', '[M-H]-', '[M+Na]+', '[M+K]+',
                                                               '[M+Cl]-', '[M+NH4]+']}
def add_adduct(s, mass_tolerance = 0.01):
    parent_mass = float(s.get("mw", None))
    s.set('parent_mass', parent_mass)
    precursor_mz = s.get("precursor_mz", None)
    for k in adducts_dict:
        ionmode_ = adducts_dict[k]['ionmode']
        charge_ = adducts_dict[k]['charge']
        correction_mass = adducts_dict[k]['correction_mass']
        if abs(precursor_mz - parent_mass - correction_mass) <= mass_tolerance:
            s = s.set('charge', charge_)
            s = s.set('ionmode', ionmode_)
            s = s.set('adduct', k)  
    return s


def apply_filters(s):
    s = add_adduct(s)
    s = default_filters(s)
    s = add_parent_mass(s, estimate_from_adduct=True)
    return s

spectrums = [apply_filters(s) for s in tqdm(spectrums) if s is not None]
np.save(os.path.join(path_data, 'preprocessed_spectrums.npy'), spectrums)


from matchms.filtering import harmonize_undefined_inchikey, harmonize_undefined_inchi, harmonize_undefined_smiles
from matchms.filtering import repair_inchi_inchikey_smiles

def clean_metadata(s):
    s = harmonize_undefined_inchikey(s)
    s = harmonize_undefined_inchi(s)
    s = harmonize_undefined_smiles(s)
    s = repair_inchi_inchikey_smiles(s)
    return s

spectrums = [clean_metadata(s) for s in tqdm(spectrums) if s is not None]


from matchms.filtering import derive_inchi_from_smiles, derive_smiles_from_inchi
from matchms.filtering import derive_inchikey_from_inchi

def clean_metadata2(s):
    s = derive_inchi_from_smiles(s)
    s = derive_smiles_from_inchi(s)
    s = derive_inchikey_from_inchi(s)
    return s

spectrums = [clean_metadata2(s) for s in tqdm(spectrums) if s is not None]
np.save(os.path.join(path_data, 'preprocessed_spectrums.npy'), spectrums)


for spectrum in tqdm(spectrums):
    name_original = spectrum.get("compound_name")
    name = name_original.replace("F dial M", "")
    # Remove last word if likely not correct:
    if name.split(" ")[-1] in ["M", "M?", "?", "M+2H/2", "MS34+Na", "M]", "Cat+M]", "Unk", "--"]:
        name = " ".join(name.split(" ")[:-1]).strip()
    if name != name_original:
        print(f"Changed compound name from {name_original} to {name}.")
        spectrum.set("compound_name", name)
        

for spec in spectrums:
    if spec.get("adduct") in ['[M+CH3COO]-/[M-CH3]-',
                             '[M-H]-/[M-Ser]-',
                             '[M-CH3]-']:
        if spec.get("ionmode") != "negative":
            spec.set("ionmode", "negative")


from matchms.filtering import normalize_intensities
from matchms.filtering import require_minimum_number_of_peaks
from matchms.filtering import select_by_mz

def post_process(s):
    s = normalize_intensities(s)
    s = select_by_mz(s, mz_from=10.0, mz_to=1000)
    s = require_minimum_number_of_peaks(s, n_required=5)
    return s

spectrums = [post_process(s) for s in tqdm(spectrums)]
spectrums = [s for s in spectrums if s is not None]


spectrums_positive = []
spectrums_negative = []
for i, spec in enumerate(spectrums):
    if spec.get("ionmode") == "positive":
        spectrums_positive.append(spec)
    elif spec.get("ionmode") == "negative":
        spectrums_negative.append(spec)
    else:
        print(f"No ionmode found for spectrum {i} ({spec.get('ionmode')})")

pickle.dump(spectrums_negative, 
            open(os.path.join(path_data, 'ALL_NIST20_negative_cleaned.pickle'), "wb"))

pickle.dump(spectrums_positive, 
            open(os.path.join(path_data, 'ALL_NIST20_positive_cleaned.pickle'), "wb"))


