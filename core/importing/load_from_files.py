# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 08:43:08 2023

@author: DELL
"""


from typing import List

from matchms.importing import load_from_msp
from matchms.importing import load_from_mgf

from core.Spectrum import Spectrum
from core.importing.load_from_mat import load_from_mat


def load_from_files(filenames: List[str]) -> List[Spectrum]:
    output = []
    for filename in filenames:
        patt = filename.split('.')[-1]
        if patt == 'msp':
            spectrums = load_from_msp(filename)
            for s in spectrums:
                output.append(Spectrum(mz=s.mz,
                                       intensities=s.intensities,
                                       metadata=s.metadata))
        elif patt == 'mgf':
            spectrums = load_from_mgf(filename)
            for s in spectrums:
                output.append(Spectrum(mz=s.mz,
                                       intensities=s.intensities,
                                       metadata=s.metadata))
        elif patt == 'mat':
            spectrums = load_from_mat(filename)
            for s in spectrums:
                output.append(s)
        else:
            continue
    return output
    

