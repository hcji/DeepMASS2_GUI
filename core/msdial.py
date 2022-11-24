# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 10:22:58 2022

@author: DELL
"""

import base64
import numpy as np
import pandas as pd
from tqdm import tqdm
from matchms import Spectrum
from core.identification import spectrum_processing


def load_MS_DIAL_Peaklist(filename, exclude_precursor = True):
    """
    Load aligned result exported by MS-DIAL and convert into a set of matchms::spectrum object.
    Arguments:
        filename: str, the path of the MS-DIAL export.
    Returns:
        List of matchms::spectrum.
    Example:
        filename = 'data/MS_DIAL_Export_CX.csv'
        load_MS_DIAL(filename)
    """
    
    if filename.split('.')[-1] == 'csv':
        data = pd.read_csv(filename)
    elif filename.split('.')[-1] == 'txt':
        data = pd.read_csv(filename, '\t')
    else:
        return None
        
    output = []
    for i in tqdm(data.index):
        s = str(data.loc[i, 'MSMS spectrum'])
        precursor_mz = float(data.loc[i, 'Precursor m/z'])
        if s == 'nan':
            continue
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
        rt = float(data.loc[i, 'RT (min)'])
        name = 'Compound_' + str(data.loc[i, 'PeakID'])
        smiles = str(data.loc[i, 'SMILES'])
        adduct = str(data.loc[i, 'Adduct'])
        isotope = str(data.loc[i, 'MS1 isotopes'])
        isotope = isotope.split(' ')
        isotope_mz = np.array([float(ss.split(':')[0]) for ss in isotope])
        isotope_intensity = np.array([float(ss.split(':')[1]) for ss in isotope])
        
        if len(mz) <= 1:
            continue
        
        obj = Spectrum(mz = mz, intensities = intensity,
                       metadata={"precursor_mz": precursor_mz,
                                 "compound_name": name,
                                 "rt": rt,
                                 "smiles": smiles,
                                 "adduct": adduct,
                                 "isotope_mz": base64.b64encode(str(isotope_mz).encode("ascii")),
                                 "isotope_intensity": base64.b64encode(str(isotope_intensity).encode("ascii"))})
        output.append(spectrum_processing(obj))
    return output



def load_MS_DIAL_Alginment(filename, exclude_precursor = True):
    """
    Load aligned result exported by MS-DIAL and convert into a set of matchms::spectrum object.
    Arguments:
        filename: str, the path of the MS-DIAL export.
    Returns:
        List of matchms::spectrum.
    Example:
        filename = 'data/MS_DIAL_Export_CX.csv'
        load_MS_DIAL(filename)
    """
    
    if filename.split('.')[-1] == 'csv':
        data = pd.read_csv(filename)
    elif filename.split('.')[-1] == 'txt':
        data = pd.read_csv(filename, '\t')
    else:
        return None
        
    output = []
    for i in tqdm(data.index):
        s = str(data.loc[i, 'MS/MS spectrum'])
        precursor_mz = float(data.loc[i, 'Average Mz'])
        if s == 'nan':
            continue
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
        name = name = 'Compound_' + str(data.loc[i, 'Alignment ID'])
        smiles = str(data.loc[i, 'SMILES'])
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
                                 "adduct": adduct,
                                 "isotope_mz": base64.b64encode(str(isotope_mz).encode("ascii")),
                                 "isotope_intensity": base64.b64encode(str(isotope_intensity).encode("ascii")),
                                 "sample_abundance": str(sample_abundance)})
        output.append(spectrum_processing(obj))
    return output