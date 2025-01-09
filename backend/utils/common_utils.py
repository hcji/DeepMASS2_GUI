import uuid

from matchms import Spectrum


def get_title_from_spectrum(spectrum: Spectrum, idx=None):
    if "compound_name" in spectrum.metadata.keys():
        return spectrum.metadata["compound_name"]
    else:
        if idx is not None:
            return f"Spectrum {idx}"
        else:
            return f"Spectrum {uuid.uuid4()}"
