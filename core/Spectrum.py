# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:53:12 2023

@author: DELL
"""


import numpy as np
from typing import Optional
from matchms.Spectrum import Fragments
from matchms.Spectrum import Spectrum as mSpectrum


class Spectrum(mSpectrum):
    
    def __init__(self, mz: np.array,
                 intensities: np.array,
                 isotopic_mz: Optional[np.array] = np.array([]),
                 isotopic_intensities: Optional[np.array] = np.array([]),
                 metadata: Optional[dict] = None,
                 metadata_harmonization: bool = True):
        
        super().__init__(mz = mz,
                         intensities = intensities,
                         metadata = metadata,
                         metadata_harmonization = metadata_harmonization)
        
        if (len(isotopic_mz) == len(isotopic_intensities)) and len(isotopic_mz) > 0:
            self.isotopic_pattern = Fragments(mz=mz, intensities=intensities)
        else:
            self.isotopic_pattern = None

