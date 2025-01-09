from matchms import Spectrum
from matchms import filtering as msfilters


def pre_process_spectrum(spectrum: Spectrum):
    spectrum = msfilters.add_parent_mass(spectrum)
    spectrum = msfilters.add_precursor_mz(spectrum)
    spectrum = msfilters.add_compound_name(spectrum)
    spectrum = msfilters.default_filters(spectrum)
    return spectrum
