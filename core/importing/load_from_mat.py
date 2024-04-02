# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 16:22:08 2023

@author: DELL
"""


from typing import Generator

import numpy as np
from matchms.importing.load_from_msp import _parse_line_with_peaks
from matchms.importing.load_from_msp import contains_metadata
from matchms.importing.load_from_msp import parse_metadata

from core.Spectrum import Spectrum


def load_from_mat(
    filename: str, metadata_harmonization: bool = True
) -> Generator[Spectrum, None, None]:
    for spectrum in parse_mat_file(filename):
        metadata = spectrum.get("params", None)
        mz = spectrum["m/z array"]
        intensities = spectrum["intensity array"]
        isotopic_mz = spectrum["isotopic m/z array"]
        isotopic_intensities = spectrum["isotopic intensity array"]
        peak_comments = spectrum["peak comments"]
        if peak_comments != {}:
            metadata["peak_comments"] = peak_comments

        # Sort by mz (if not sorted already)
        if not np.all(mz[:-1] <= mz[1:]):
            idx_sorted = np.argsort(mz)
            mz = mz[idx_sorted]
            intensities = intensities[idx_sorted]

        yield Spectrum(
            mz=mz,
            intensities=intensities,
            isotopic_mz=isotopic_mz,
            isotopic_intensities=isotopic_intensities,
            metadata=metadata,
            metadata_harmonization=metadata_harmonization,
        )


def parse_mat_file(filename: str) -> Generator[dict, None, None]:
    """Read msp file and parse info in List of spectrum dictionaries."""

    # Lists/dicts that will contain all params, masses and intensities of each molecule
    params = {}
    masses = np.array([])
    intensities = np.array([])
    isotopic_masses = np.array([])
    isotopic_intensities = np.array([])
    peak_comments = {}
    params["mstype"] = "MS1"

    # Peaks counter. Used to track and count the number of peaks
    peakscount = 0

    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            rline = line.rstrip()

            if len(rline) == 0:
                continue

            if contains_metadata(rline):
                parse_metadata(rline, params)
                continue

            mz, ints, comment = _parse_line_with_peaks(rline)

            if params["mstype"] == "MS2":
                masses = np.append(masses, mz)
                intensities = np.append(intensities, ints)
                peakscount += len(mz)
            elif params["mstype"] == "MS1":
                isotopic_masses = np.append(isotopic_masses, mz)
                isotopic_intensities = np.append(isotopic_intensities, ints)
            else:
                pass

            if comment is not None:
                peak_comments.update({masses[-1]: comment})

            # Obtaining the masses and intensities
            if (int(params["num peaks"]) == peakscount) and (params["mstype"] == "MS2"):
                peakscount = 0
                yield {
                    "params": (params),
                    "m/z array": masses,
                    "intensity array": intensities,
                    "isotopic m/z array": isotopic_masses,
                    "isotopic intensity array": isotopic_intensities,
                    "peak comments": peak_comments,
                }

                params = {}
                masses = []
                intensities = []
                isotopic_masses = []
                isotopic_intensities = []
                peak_comments = {}
                params["mstype"] = "MS1"
