# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:53:12 2023

@author: DELL
"""


from typing import Optional

import numpy as np
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
            self.isotopic_pattern = Fragments(mz=isotopic_mz, intensities=isotopic_intensities)
            self.clean_isotopic_pattern()
        else:
            self.isotopic_pattern = None


    def clean_isotopic_pattern(self, tolerence = 0.01):
        precursor_mz = self.get('precursor_mz')
        isotopic_mz_ = precursor_mz + np.arange(5) * 1.00866
        minimum_diff = np.array([np.min(np.abs(m - isotopic_mz_)) for m in self.isotopic_pattern.mz])
        k = np.where(minimum_diff <= tolerence)[0]
        if len(k) >= 2:
            self.isotopic_pattern = Fragments(mz=self.isotopic_pattern.mz, intensities=self.isotopic_pattern.intensities)
        else:
            self.isotopic_pattern = None
