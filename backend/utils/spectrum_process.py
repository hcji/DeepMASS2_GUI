from matchms import Spectrum
from matchms.filtering import (
    default_filters,
    normalize_intensities,
    add_parent_mass,
    correct_charge,
)
from matchms.importing import (
    load_from_msp,
    load_from_mgf,
    load_from_mzml,
    load_from_mzxml,
    load_from_json,
)


def load_spectrum_file(file_name: str):
    '''
    从文件名中读取质谱
    Args:
        file_name:

    Returns:

    '''
    if file_name.lower().endswith("msp"):
        spectra_list = load_from_msp(file_name)
    elif file_name.lower().endswith("mgf"):
        spectra_list = load_from_mgf(file_name)
    elif file_name.lower().endswith("mzml"):
        spectra_list = load_from_mzml(file_name)
    elif file_name.lower().endswith("mzxml"):
        spectra_list = load_from_mzxml(file_name)
    elif file_name.lower().endswith("json"):
        spectra_list = load_from_json(file_name)
    # TODO 写mat文件读取方法
    elif file_name.lower().endswith("mat"):
        spectra_list = []
    else:
        spectra_list = []
    spectra_list = [spectra for spectra in spectra_list]
    if len(spectra_list) > 0:
        spectra_list = [deepmass_default_filter(spectrum) for spectrum in spectra_list]
    return spectra_list


def deepmass_default_filter(spectrum: Spectrum):
    spectrum = default_filters(spectrum)
    spectrum = correct_charge(spectrum)
    spectrum = normalize_intensities(spectrum)
    spectrum = add_parent_mass(spectrum)
    return spectrum