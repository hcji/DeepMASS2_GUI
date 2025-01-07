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


def load_spectrum_file(filename: str):
    """
    从文件名中读取质谱
    Args:
        file_name:

    Returns:

    # """
    # if file_name.lower().endswith("msp"):

    #     spectra_list = load_from_msp(file_name)

    # elif file_name.lower().endswith("mgf"):
    #     spectra_list = load_from_mgf(file_name)

    # # 写mat文件读取方法
    # elif file_name.lower().endswith("mat"):
    #     spectra_list = load_from_mat(file_name)
    # else:
    #     spectra_list = []
    # spectra_list = [spectra for spectra in spectra_list]
    # if len(spectra_list) > 0:
    #     spectra_list = [deepmass_default_filter(spectrum) for spectrum in spectra_list]
    # return spectra_list

    output = []
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
        pass
    output_1 = []
    for s in output:
        if s.get('compound_name') is None:
            s.set('compound_name', 'Unknown_{}'.format(len(output_1)))
        output_1.append(s)
    return output_1


def fill_pos_charge(spectrum: Spectrum):
    """
    当没有识别到模式时，默认注入正模式
    Args:
        spectrum:

    Returns:

    """
    if (
        spectrum.metadata.get("ionmode") == "n/a"
        or spectrum.metadata.get("ionmode") == None
    ):
        spectrum.set("ionmode", "positive")
    return spectrum

if __name__ == "__main__":
    print(load_spectrum_file("./example/minimum_example(2).mat"))
