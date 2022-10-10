# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 14:55:39 2022

@author: DELL
"""


import numpy as np
import pandas as pd
from tqdm import tqdm

from matchms.Spectrum import Spectrum
from matchms.similarity import CosineGreedy

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS


def disable_rdkit_logging():
    """
    Disables RDKit whiny logging.
    """
    import rdkit.rdBase as rkrb
    import rdkit.RDLogger as rkl
    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog('rdApp.error')

disable_rdkit_logging()
get_fp = lambda x: AllChem.GetMorganFingerprintAsBitVect(x, radius=2)
get_sim = lambda x, y: DataStructs.DiceSimilarity(x, y)


def load_MS_DIAL(filename, exclude_reduplicated = True, exclude_precursor = False):
    """
    Load aligned result exported by MS-DIAL and convert into a set of matchms::spectrum object.
    Arguments:
        filename: str, the path of the MS-DIAL export.
    Returns:
        List of matchms::spectrum.
    Example:
        filename = 'data/MS_DIAL_Export_MTBLS4099_POS.csv'
        load_MS_DIAL(filename)
    """
    data = pd.read_csv(filename)
    output, all_smiles = [], []
    for i in tqdm(data.index):
        s = str(data.loc[i, 'MS/MS spectrum'])
        precursor_mz = float(data.loc[i, 'Average Mz'])
        if s == 'nan':
            mz = np.array([])
            intensity = np.array([])
        else:
            s = s.split(' ')
            mz = np.array([float(ss.split(':')[0]) for ss in s if ':' in ss])
            intensity = np.array([float(ss.split(':')[1]) for ss in s if ':' in ss])      
            if exclude_precursor:
                k = np.where(np.logical_and(mz <= precursor_mz - 0.1, intensity > 0))[0]
            else:
                k = np.where(intensity > 0)[0]
            mz = mz[k]
            intensity = intensity[k]
            intensity /= (np.max(intensity) + 10 **-10)
        sample_cols = [s for s in data.columns if 'Sample ' in s]
        rt = float(data.loc[i, 'Average Rt(min)'])
        name = str(data.loc[i, 'Metabolite name'])
        smiles = str(data.loc[i, 'SMILES'])
        avg_abund = float(data.loc[i, 'Average'])
        if smiles != 'nan':
            m = Chem.MolFromSmiles(smiles)
            if m is not None:
                smiles = Chem.MolToSmiles(m)
            else:
                smiles = 'nan'
                
        exact_mass = precursor_mz - 1.0000
        adduct = str(data.loc[i, 'Adduct type'])
        isotope = str(data.loc[i, 'MS1 isotopic spectrum'])
        isotope = isotope.split(' ')
        isotope_mz = np.array([float(ss.split(':')[0]) for ss in isotope])
        isotope_intensity = np.array([float(ss.split(':')[1]) for ss in isotope])
        sample_abundance = np.array(data.loc[i, sample_cols])
        
        obj = Spectrum(mz = mz, intensities = intensity,
                       metadata={"precursor_mz": precursor_mz,
                                 "compound_name": name,
                                 "rt": rt,
                                 "smiles": smiles,
                                 "adduct_type": adduct,
                                 "isotope_mz": isotope_mz,
                                 "isotope_intensity": isotope_intensity,
                                 "sample_abundance": sample_abundance,
                                 "average_abundance": avg_abund})
        
        if exclude_reduplicated:
            if smiles in all_smiles:
                wh = all_smiles.index(smiles)
                if output[wh][-1].metadata['average_abundance'] >= avg_abund:
                    continue
                else:
                    output[wh][-1] = obj
                    continue
         
        if smiles != 'nan':
            all_smiles.append(smiles)
            score = 1
            reference = 'Explicit'
        else:
            score = 0
            reference = ''
        
        output.append([precursor_mz, adduct, exact_mass, smiles, score, reference, obj])
    output = pd.DataFrame(output, columns=['precursor_mz', 'adduct', 'exact_mass', 'smiles', 'score', 'reference', 'spectrum'])
    return output


def get_tagged_atoms_from_mol(mol):
    '''Takes an RDKit molecule and returns list of tagged atoms and their
    corresponding numbers'''
    atoms = []
    atom_tags = []
    for atom in mol.GetAtoms():
        if atom.HasProp('molAtomMapNumber'):
            atoms.append(atom)
            atom_tags.append(int(atom.GetProp('molAtomMapNumber')))
    return atom_tags


def calc_frag_mass(frag):
    '''Takes an RDKit fragment and returns the exact mass of the fragment'''
    mass = 0
    for ad in frag.GetAtoms():
        mass += ad.GetMass()
    return mass
    

def calc_possible_spectrum_loss(smiles_1, smiles_2):
    """
    Calculate mass difference between related fragments of two compounds.
    Arguments:
        smiles_1, smiles_2: str, two different smiles of compounds.
    Returns:
        DataFrame, 
            transform, transformation from smiles_1 to smiles_2.
            loss, corresponding neutral loss of the transformations.
    Example:
        smiles_1 = 'COc1cc(O)c2c(c1)OC(c1ccc(O)cc1)CC2=O'
        smiles_2 = 'CC1OC(OCC2OC(Oc3cc(O)c4c(c3)OC(c3ccc(O)cc3)CC4=O)C(O)C(O)C2O)C(O)C(O)C1O'
        calc_possible_spectrum_loss(smiles_1, smiles_2)
    """    
    
    try:
        x = Chem.MolToSmiles(Chem.MolFromSmiles(smiles_1))
        y = Chem.MolToSmiles(Chem.MolFromSmiles(smiles_2))
    except:
        return None
    
    mol1 = Chem.AddHs(Chem.MolFromSmiles(x))
    mol2 = Chem.AddHs(Chem.MolFromSmiles(y))
    
    if get_sim(get_fp(mol1), get_fp(mol2)) < 0.3:
        return None
    
    mcs = rdFMCS.FindMCS([mol1, mol2], bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                         matchValences = True, ringMatchesRingOnly = True)
    if mcs.numAtoms <= 5:
        return None
    
    mcs_str = mcs.smartsString
    
    rdu1 = AllChem.DeleteSubstructs(mol1, Chem.MolFromSmarts(mcs_str))
    rdu2 = AllChem.DeleteSubstructs(mol2, Chem.MolFromSmarts(mcs_str))
    
    try:
        rdu1 = Chem.GetMolFrags(rdu1, asMols=True)
    except:
        rdu1 = np.array([rdu1])
        
    try:
        rdu2 = Chem.GetMolFrags(rdu2, asMols=True)
    except:
        rdu2 = np.array([rdu2])
    
    if (len(rdu1) == 0) and (len(rdu2) == 0):
        return None

    mass_1 = np.array([calc_frag_mass(m) for m in rdu1])
    mass_2 = np.array([calc_frag_mass(m) for m in rdu2])
    if len(mass_1) == 0:
        mass_1 = np.array([0])
    if len(mass_2) == 0:
        mass_2 = np.array([0])
    
    mass_diffs, mol_transform = [], []
    for i in range(len(mass_1)):
        for j in range(len(mass_2)):
            try:
                a = Chem.MolToSmiles(Chem.RemoveHs(rdu1[i]))
            except:
                a = 'None'
            try:
                b = Chem.MolToSmiles(Chem.RemoveHs(rdu2[j]))
            except:
                b = 'None'
            mol_transform.append('{}>>{}'.format(a,b))
            mass_diffs.append(mass_2[j] - mass_1[i])
    
    return pd.DataFrame({'transform': mol_transform, 'loss': mass_diffs})
    
    
def calc_aligned_similarity(smiles_1, smiles_2, spectrum_1, spectrum_2, mz_tol=0.05, similarity_function=CosineGreedy()):
    """
    Calculate dtw similarity between two spectrums.
    Arguments:
        smiles_1, smiles_2: str, two different smiles of compounds.
        spectrum_1, spectrum_2: Two different spectrum of matchms.
    Returns:
        similarity: float, similarity between aligned spectrums.
        matching_data: DataFrame, fragment matching information.
    Example:
        smiles_1 = 'CCCC=C1C2=CC=CC=C2C(=O)O1'
        smiles_2 = 'CCCC=C1C2=C(C=CCC2)C(=O)O1'
        spectrum_1 = Spectrum(mz = np.array([91.1, 115.1, 117.1, 128.1, 129.1, 143.1, 145.1, 152.1, 153.1, 171.1, 189.1]), 
                              intensities = np.array([0.12314933, 0.10446688, 0.16478671, 0.56083889, 0.11087135,
                                                      0.43528005, 0.1149675 , 0.10339803, 0.51058281, 0.999999, 0.88490263]),
                              metadata={"precursor_mz": 189.0909})
        spectrum_2 = Spectrum(mz = np.array([ 79.1, 93.1, 105.1, 117.1, 145.1, 173.1, 191.1]), 
                              intensities = np.array([0.10704697, 0.10657389, 0.1382483 , 0.12679477, 0.16397634,
                                                      0.26150501, 0.999999]),
                              metadata={"precursor_mz": 191.1064})
        calc_aligned_similarity(smiles_1, smiles_2, spectrum_1, spectrum_2)
    """
    
    loss = calc_possible_spectrum_loss(smiles_1, smiles_2)
    
    if loss is None:
        loss_1 = loss_2 = 0
        loss = pd.DataFrame({'transform': [], 'loss': []})
    else:
        loss_1 = -sum([l for l in list(loss['loss']) if l < 0])
        loss_2 = sum([l for l in list(loss['loss']) if l > 0])
    
    mcs1 = 9999
    mcs2 = 9999
    try:
        mcs1 = spectrum_1.metadata['precursor_mz'] - loss_1
    except:
        pass
    try:
        mcs2 = spectrum_2.metadata['precursor_mz'] - loss_2
    except:
        pass
    maxCS = min(mcs1, mcs2)
    
    if (len(spectrum_1.mz) == 0) or (len(spectrum_2.mz) == 0):
        return 0, None
    
    x_mz, x_intensities = spectrum_1.mz, spectrum_1.intensities
    y_mz, y_intensities = spectrum_2.mz, spectrum_2.intensities
    
    y_mz_new, y_intensities_new = [], []
    matching_data = []
    for i, y_mz_ in enumerate(y_mz):
        if y_intensities[i] < 0.01:
            continue
        if np.min(np.abs(y_mz_ - x_mz)) <= mz_tol:
            if y_mz_ > maxCS + 2.006:
                continue
            a = y_mz_
            b = x_mz[np.argmin(np.abs(y_mz_ - x_mz))]
            c = abs(a - b)
            d = y_intensities[i]
            y_mz_new.append(y_mz_)
            y_intensities_new.append(y_intensities[i])
            matching_data.append([a, b, c, d])
        else:
            matched = False
            for loss_ in loss['loss']:
                '''
                if y_mz_ - loss_ > maxCS - loss_ + 2.006:
                    continue
                '''
                if np.min(np.abs(y_mz_ - loss_ - x_mz)) <= mz_tol:
                    matched = True
                    a = y_mz_
                    b = x_mz[np.argmin(np.abs(y_mz_ - loss_ - x_mz))]
                    c = abs(a - b)
                    d = y_intensities[i]
                    y_mz_new.append(y_mz_ - loss_)
                    y_intensities_new.append(y_intensities[i])
                    matching_data.append([a, b, c, d])
                    break
            if not matched:
                y_mz_new.append(y_mz_)
                y_intensities_new.append(y_intensities[i])
    y_mz_new = np.array(y_mz_new)
    y_intensities_new = np.array(y_intensities_new)
    
    index = np.argsort(y_mz_new)
    y_mz_new = y_mz_new[index]
    y_intensities_new = y_intensities_new[index]
    
    spectrum_2_aligned = Spectrum(mz = y_mz_new, 
                                  intensities = y_intensities_new,
                                  metadata = spectrum_2.metadata)
    
    similarity = float(similarity_function.pair(spectrum_1, spectrum_2_aligned)['score'])
    matching_data = pd.DataFrame(matching_data, columns = ['reference', 'query', 'loss', 'intensity'])
    return similarity, matching_data
