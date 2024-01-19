# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 10:57:51 2023

@author: DELL
"""

import base64

import numpy as np
import pandas as pd
from matchms import Spectrum
from tqdm import tqdm

from core.identification import spectrum_processing


def load_MS_DIAL_Peaklist(filename, exclude_precursor=False):
    """
    Load aligned result exported by MS-DIAL and convert into a set of matchms::spectrum object.
    Arguments:
        filename: str, the path of the MS-DIAL export.
    Returns:
        List of matchms::spectrum.
    """

    if filename.split(".")[-1] == "csv":
        data = pd.read_csv(filename)
    elif filename.split(".")[-1] == "txt":
        data = pd.read_csv(filename, "\t")
    else:
        return None

    output = []
    for i in tqdm(data.index):
        s = str(data.loc[i, "MSMS spectrum"])
        precursor_mz = float(data.loc[i, "Precursor m/z"])
        if s == "nan":
            continue
        else:
            s = s.split(" ")
            mz = np.array([float(ss.split(":")[0]) for ss in s if ":" in ss])
            intensity = np.array([float(ss.split(":")[1]) for ss in s if ":" in ss])
            if exclude_precursor:
                k = np.where(np.logical_and(mz <= precursor_mz - 0.1, intensity > 0))[0]
            else:
                k = np.where(intensity > 0)[0]
            mz = mz[k]
            intensity = intensity[k]
            intensity /= np.max(intensity) + 10**-10

            ww = np.where(intensity >= 0.05)[0]
            mz = mz[ww]
            intensity = intensity[ww]

        rt = float(data.loc[i, "RT (min)"])
        index = "Peak_" + str(data.loc[i, "PeakID"])
        smiles = str(data.loc[i, "SMILES"])
        adduct = str(data.loc[i, "Adduct"])
        isotope = str(data.loc[i, "MS1 isotopes"])
        isotope = isotope.split(" ")
        isotope_mz = np.array([float(ss.split(":")[0]) for ss in isotope])
        isotope_intensity = np.array([float(ss.split(":")[1]) for ss in isotope])

        if len(mz) <= 1:
            continue

        obj = Spectrum(
            mz=mz,
            intensities=intensity,
            metadata={
                "precursor_mz": precursor_mz,
                "peak_index": index,
                "rt": rt,
                "smiles": smiles,
                "adduct": adduct,
                "isotope_mz": base64.b64encode(str(isotope_mz).encode("ascii")),
                "isotope_intensity": base64.b64encode(
                    str(isotope_intensity).encode("ascii")
                ),
            },
        )
        output.append(spectrum_processing(obj))
    return output


def load_MS_DIAL_Alginment(filename, exclude_precursor=False, sample_cols=[]):
    """
    Load aligned result exported by MS-DIAL and convert into a set of matchms::spectrum object.
    Arguments:
        filename: str, the path of the MS-DIAL export.
    Returns:
        List of matchms::spectrum.
    Example:
        filename = 'example/Plasma/ms_dial_positive.csv'
        load_MS_DIAL_Alginment(filename)
    """

    if filename.split(".")[-1] == "csv":
        data = pd.read_csv(filename)
    elif filename.split(".")[-1] == "txt":
        data = pd.read_csv(filename, "\t")
    else:
        return None

    output = []
    for i in tqdm(data.index):
        s = str(data.loc[i, "MS/MS spectrum"])
        precursor_mz = float(data.loc[i, "Average Mz"])
        if s == "nan":
            continue
        else:
            s = s.split(" ")
            mz = np.array([float(ss.split(":")[0]) for ss in s if ":" in ss])
            intensity = np.array([float(ss.split(":")[1]) for ss in s if ":" in ss])
            if exclude_precursor:
                k = np.where(np.logical_and(mz <= precursor_mz - 0.1, intensity > 0))[0]
            else:
                k = np.where(intensity > 0)[0]
            mz = mz[k]
            intensity = intensity[k]
            intensity /= np.max(intensity) + 10**-10

            ww = np.where(intensity >= 0.05)[0]
            mz = mz[ww]
            intensity = intensity[ww]

        rt = float(data.loc[i, "Average Rt(min)"])
        index = "Peak_" + str(data.loc[i, "Alignment ID"])
        smiles = str(data.loc[i, "SMILES"])
        adduct = str(data.loc[i, "Adduct type"])
        isotope = str(data.loc[i, "MS1 isotopic spectrum"])
        isotope = isotope.split(" ")
        isotope_mz = np.array([float(ss.split(":")[0]) for ss in isotope])
        isotope_intensity = np.array([float(ss.split(":")[1]) for ss in isotope])
        sample_abundance = np.array(data.loc[i, sample_cols])
        precursor_intensity = np.nanmean(sample_abundance)

        obj = Spectrum(
            mz=mz,
            intensities=intensity,
            metadata={
                "precursor_mz": precursor_mz,
                "peak_index": index,
                "rt": rt,
                "smiles": smiles,
                "adduct": adduct,
                "precursor_intensity": precursor_intensity,
                "isotope_mz": base64.b64encode(str(isotope_mz).encode("ascii")),
                "isotope_intensity": base64.b64encode(
                    str(isotope_intensity).encode("ascii")
                ),
            },
        )
        output.append(spectrum_processing(obj))
    return output


def remove_duplicate(spectrums):
    new_spectrums = []
    rt, mz, iontype, intensities = [], [], [], []
    for s in tqdm(spectrums):
        [rt_, mz_, iontype_, intensity_, adduct_] = [
            s.metadata[k]
            for k in [
                "retention_time",
                "precursor_mz",
                "ionmode",
                "precursor_intensity",
                "adduct",
            ]
        ]
        if adduct_ not in ["[M+H]+", "[M-H]-"]:
            continue
        wh = np.logical_and(
            np.abs(np.array(rt) - rt_) < 18,
            np.abs(np.array(mz) - mz_) < 0.01,
            np.array([i == iontype_ for i in iontype]),
        )
        wh = np.where(wh)[0]
        if len(wh) > 0:
            w = wh[0]
            if intensity_ >= intensities[w]:
                new_spectrums[w] = s
                intensities[w] = intensity_
            else:
                continue
        else:
            rt.append(rt_)
            mz.append(mz_)
            iontype.append(iontype_)
            intensities.append(intensity_)
            new_spectrums.append(spectrum_processing(s))
    return new_spectrums


def save_as_sirius(spectrums, export_path):
    for i, s in enumerate(spectrums):
        energy = 35
        compound = s.get("compound_name")
        parentmass = s.get("parent_mass")
        ionization = s.get("adduct")

        isotope_mz = (
            base64.b64decode(s.get("isotope_mz")).decode("ascii").replace("\n", "")
        )
        isotope_intensity = (
            base64.b64decode(s.get("isotope_intensity"))
            .decode("ascii")
            .replace("\n", "")
        )
        isotope_mz = [
            float(s)
            for s in isotope_mz.replace("[", "").replace("]", "").split(" ")
            if s != ""
        ]
        isotope_intensity = [
            float(s)
            for s in isotope_intensity.replace("[", "").replace("]", "").split(" ")
            if s != ""
        ]

        with open(export_path + "/" + compound + ".ms", "w") as ms:
            ms.write(">compound {}\n".format(compound))
            ms.write(">ionization {}\n".format(ionization))
            ms.write("\n")

            ms.write(">collision {}\n".format(energy))

            for p in range(len(s.mz)):
                mz = s.mz[p]
                intensity = s.intensities[p]
                ms.write("{} {}\n".format(mz, intensity))

            ms.write("\n\n")
            ms.write(">ms1peaks\n")
            for p in range(len(isotope_mz)):
                mz = isotope_mz[p]
                intensity = isotope_intensity[p]
                ms.write("{} {}\n".format(mz, intensity))
            ms.write("\n")
    pass


def save_as_msfinder(spectrums, export_path):
    for i, s in enumerate(spectrums):
        compound = s.get("compound_name")
        precursor_mz = s.get("precursor_mz")
        ionmode = s.get("ionmode").capitalize()
        ionization = s.get("adduct")

        isotope_mz = (
            base64.b64decode(s.get("isotope_mz")).decode("ascii").replace("\n", "")
        )
        isotope_intensity = (
            base64.b64decode(s.get("isotope_intensity"))
            .decode("ascii")
            .replace("\n", "")
        )
        isotope_mz = [
            float(s)
            for s in isotope_mz.replace("[", "").replace("]", "").split(" ")
            if s != ""
        ]
        isotope_intensity = [
            float(s)
            for s in isotope_intensity.replace("[", "").replace("]", "").split(" ")
            if s != ""
        ]

        with open(export_path + "/" + compound + ".mat", "w") as ms:
            ms.write("NAME: {}\n".format(compound))
            ms.write("PRECURSORMZ: {}\n".format(precursor_mz))
            ms.write("PRECURSORTYPE: {}\n".format(ionization))
            ms.write("IONMODE: {}\n".format(ionmode))
            ms.write("\n")

            ms.write("MSTYPE: MS1\n")
            ms.write("Num Peaks: {}\n".format(len(isotope_mz)))
            for p in range(len(isotope_mz)):
                mz = isotope_mz[p]
                intensity = isotope_intensity[p]
                ms.write("{}\t{}\n".format(mz, intensity))
            ms.write("\n")

            ms.write("MSTYPE: MS2\n")
            ms.write("Num Peaks: {}\n".format(len(s.mz)))
            for p in range(len(s.mz)):
                mz = s.mz[p]
                intensity = s.intensities[p]
                ms.write("{}\t{}\n".format(mz, intensity))

            ms.write("\n")
    pass


def spectrum_to_vector(s, min_mz=0, max_mz=1000, scale=0.1):
    """
    Convert spectrum object to vector.
    Arguments:
        s: matchms::spectrum
        min_mz: float, start of mz value.
        max_mz: float, end of mz value.
        scale: float, scale of mz bin.
    Returns:
        Numpy array of spectrum.
    """
    bit = round((1 + max_mz - min_mz) / scale)
    vec = np.zeros(bit)
    if s is None:
        return vec
    else:
        k = np.logical_and(min_mz <= s.mz, s.mz <= max_mz)
        idx = np.round((s.mz[k] - min_mz) / scale).astype(int)
        val = s.intensities[k]
        vec[idx] = val
        vec = vec / (np.max(vec) + 10**-6)
        return vec


def consensus_spectrum(spectrums, mz_window=0.2):
    tot_array = []
    for i, s in enumerate(spectrums):
        mz, intensity = s.peaks.mz, s.peaks.intensities
        array = np.vstack((mz, intensity, np.repeat(i, len(mz)))).T
        tot_array.append(array)

    i = 0
    mz, intensity = [], []
    tot_array = np.vstack(tot_array)
    while True:
        if i >= len(tot_array):
            break
        m = tot_array[i, 0]
        j = np.searchsorted(tot_array[:, 0], m + mz_window)
        a = tot_array[i:j, 0]
        b = tot_array[i:j, 1]
        a = np.round(np.sum(a * b) / np.sum(b), 5)
        b = np.round(np.max(b), 5)
        mz.append(a)
        intensity.append(b)
        i = j
    output = np.vstack((mz, intensity)).T
    return output
