# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 09:34:04 2022

@author: DELL
"""


from matchms.exporting import save_as_mgf
from matchms.importing import load_from_mgf
from core.identification import spectrum_processing
from core.msdial import load_MS_DIAL_Peaklist

'''
mix_msdial = load_MS_DIAL_Peaklist('example/600mix_pos.txt')
mix_msdial = [spectrum_processing(s) for s in mix_msdial]
mix_msdial = [s for s in mix_msdial if s is not None]
save_as_mgf(mix_msdial, 'example/600mix_pos.mgf')
'''

spectrums = [s for s in load_from_mgf('example/600mix_pos.mgf')]
