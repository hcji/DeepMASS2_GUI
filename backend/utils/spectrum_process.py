from matchms import Spectrum
from matchms.filtering import (
    default_filters,
    normalize_intensities,
    add_parent_mass,
    correct_charge,
    derive_ionmode,
)
from matchms.importing import (
    load_from_msp,
    load_from_mgf,
    load_from_mzml,
    load_from_mzxml,
    load_from_json,
)

from core.importing.load_from_mat import load_from_mat


def load_spectrum_file(file_name: str):
    """
    从文件名中读取质谱
    Args:
        file_name:

    Returns:

    """
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
    # 写mat文件读取方法
    elif file_name.lower().endswith("mat"):
        spectra_list = load_from_mat(file_name)
    else:
        spectra_list = []
    spectra_list = [spectra for spectra in spectra_list]
    if len(spectra_list) > 0:
        spectra_list = [deepmass_default_filter(spectrum) for spectrum in spectra_list]
    return spectra_list


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


def deepmass_default_filter(spectrum: Spectrum):
    spectrum = default_filters(spectrum)
    spectrum = normalize_intensities(spectrum)
    spectrum = derive_ionmode(spectrum)
    spectrum = correct_charge(spectrum)
    spectrum = fill_pos_charge(spectrum)
    spectrum = add_parent_mass(spectrum)
    return spectrum


if __name__ == "__main__":
    print(load_spectrum_file("./example/minimum_example(2).mat"))
