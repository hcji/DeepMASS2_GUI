# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 08:43:08 2023

@author: DELL
"""


from typing import List

import matchms.filtering as msfilters
from matchms.importing import load_from_msp
from matchms.importing import load_from_mgf

from core.Spectrum import Spectrum
from core.importing.load_from_mat import load_from_mat



def clean_spectrum(s):
    s = msfilters.default_filters(s)
    s = msfilters.add_compound_name(s)
    s = msfilters.correct_charge(s)
    s = msfilters.add_parent_mass(s)
    s = msfilters.normalize_intensities(s)
    return s


def load_from_files(filenames: List[str]) -> List[Spectrum]:
    output = []
    for filename in filenames:
        patt = filename.split('.')[-1]
        if patt == 'msp':
            spectrums = load_from_msp(filename)
            for s in spectrums:
                s = clean_spectrum(s)
                output.append(Spectrum(mz=s.mz,
                                       intensities=s.intensities,
                                       metadata=s.metadata))
        elif patt == 'mgf':
            spectrums = load_from_mgf(filename)
            for s in spectrums:
                s = clean_spectrum(s)
                output.append(Spectrum(mz=s.mz,
                                       intensities=s.intensities,
                                       metadata=s.metadata))
        elif patt == 'mat':
            spectrums = load_from_mat(filename)
            for s in spectrums:
                isotopic_mz = s.isotopic_pattern.mz
                isotopic_intensities = s.isotopic_pattern.intensities
                s = clean_spectrum(s)
                output.append(Spectrum(mz=s.mz,
                                       intensities=s.intensities,
                                       isotopic_mz=isotopic_mz,
                                       isotopic_intensities=isotopic_intensities,
                                       metadata=s.metadata))
        else:
            continue
    output_1 = []
    for s in output:
        if s.get('compound_name') is None:
            s.set('compound_name', 'Unknown_{}'.format(len(output_1)))
        output_1.append(s)
    return output_1
    

