# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 09:00:33 2024

@author: DELL
"""


from matchms.exporting import save_as_mgf
from core.external import load_MS_DIAL_Alginment

spectrums = load_MS_DIAL_Alginment('example/MS_DIAL_example.txt')
save_as_mgf(spectrums, 'example/MS_DIAL_example.mgf')