import heapq
import logging
from itertools import chain

import gradio as gr
import numpy as np
import pandas as pd
from matchms import Spectrum
from matchms import filtering as msfilters
from matplotlib.figure import Figure
from molmass import Formula
from rdkit import Chem
from rdkit.Chem import rdFMCS, Draw

from backend.load_config import GLOBAL_CONFIG

# 从配置文件获取质谱图DPI和比例
dpi_config = GLOBAL_CONFIG["identification"]["plot"]["dpi"]
width_config = GLOBAL_CONFIG["identification"]["plot"]["width"]
length_config = GLOBAL_CONFIG["identification"]["plot"]["length"]


_annotation_kws = {
    "horizontalalignment": "left",  # if not mirror_intensity else "right",
    "verticalalignment": "center",
    "fontsize": 2,
    "rotation": 90,
    "rotation_mode": "anchor",
    "zorder": 5,
}
_annotation_reverse_kws = {
    "horizontalalignment": "left",  # if not mirror_intensity else "right",
    "verticalalignment": "center",
    "fontsize": 2,
    "rotation": 270,
    "rotation_mode": "anchor",
    "zorder": 5,
}


def show_mol(structure_state, cur_spectrum, evt: gr.SelectData):
    # ① 空值保护
    if cur_spectrum is None or structure_state is None:
        return None, None
    
    try:
        if evt is None or evt.index is None:
            return None, None

        line_num = evt.index[0]
        return show_default_mol(structure_state, cur_spectrum, line_num)
    except Exception as e:
        logging.error(f"Error in show_mol: {e}")
        return None, None

def show_default_mol(structure_state, cur_spectrum, idx=0):
    # ① 空值保护
    if cur_spectrum is None or not hasattr(cur_spectrum, "metadata"):
        return None, None
    
    try:
        ref_smi = cur_spectrum.metadata["reference"][idx].metadata["smiles"]
        anno_img, ref_img = plot_2_mol(structure_state, ref_smi)
        return anno_img, ref_img
    except Exception as e:
        logging.error(f"Error in show_default_mol: {e}")
        return None, None

def get_formula_mass(formula: str):
    f = Formula(formula.replace("-", ""))
    try:
        mass = f.isotope.mass
    except Exception as e:
        logging.warn(f"{e}")
        mass = 0
    return mass


def add_topk_mz_text(ax, mz, intensities, is_reverse=False, top_k=3):
    nlargest_index = heapq.nlargest(
        top_k, range(len(intensities)), intensities.__getitem__
    )
    print(nlargest_index)
    for idx in nlargest_index:
        if not is_reverse:
            ax.text(mz[idx], intensities[idx], f"{round(mz[idx],4)}", _annotation_kws)
        else:
            ax.text(
                mz[idx],
                -1 * intensities[idx],
                f"{round(mz[idx],4)}",
                _annotation_reverse_kws,
            )


def plot_2_spectrum(spectrum: Spectrum, reference: Spectrum, loss=False):
    mz, abundance = spectrum.peaks.mz, spectrum.peaks.intensities
    mz1, abunds1 = reference.peaks.mz, reference.peaks.intensities
    if loss:
        try:
            spectrum = msfilters.add_parent_mass(spectrum)
            spectrum = msfilters.add_losses(
                spectrum, loss_mz_from=10.0, loss_mz_to=2000.0
            )
            reference = msfilters.add_parent_mass(reference)
            reference = msfilters.add_losses(
                reference, loss_mz_from=10.0, loss_mz_to=2000.0
            )
            mz, abundance = spectrum.losses.mz, spectrum.losses.intensities
            mz1, abunds1 = reference.losses.mz, reference.losses.intensities
        except:
            print("Cannot Plot Losses")
            abundance /= np.max(abundance)
    abunds1 /= np.max(abunds1)

    fig = Figure(figsize=(width_config, length_config), dpi=dpi_config)
    fig.subplots_adjust(top=0.95, bottom=0.3, left=0.18, right=0.95)

    axes = fig.add_subplot(111)
    axes.tick_params(width=0.8, labelsize=3)
    axes.spines["bottom"].set_linewidth(0.5)
    axes.spines["left"].set_linewidth(0.5)
    axes.spines["right"].set_linewidth(0.5)
    axes.spines["top"].set_linewidth(0.5)
    axes.tick_params(width=0.8, labelsize=3)
    axes.vlines(mz, ymin=0, ymax=abundance, color="r", lw=0.5)
    axes.vlines(mz1, ymin=0, ymax=-abunds1, color="b", lw=0.5)
    axes.axhline(y=0, color="black", lw=0.5)
    axes.set_xlabel("m/z", fontsize=3.5)
    axes.set_ylabel("abundance", fontsize=3.5)
    add_topk_mz_text(axes, mz, abundance)
    add_topk_mz_text(axes, mz1, abunds1, is_reverse=True)
    axes.set_ylim(-1.9, 1.9)
    return fig


def plot_2_mol(smi_anno, smi_ref, hightlight=True):
    mol_anno = Chem.MolFromSmiles(smi_anno)
    mol_ref = Chem.MolFromSmiles(smi_ref)
    if hightlight:
        mcs = rdFMCS.FindMCS(
            [mol_anno, mol_ref],
            bondCompare=rdFMCS.BondCompare.CompareOrderExact,
            matchValences=True,
            ringMatchesRingOnly=True,
        )
        mcs_str = mcs.smartsString
        mcs_mol = Chem.MolFromSmarts(mcs_str)
        all_subs_anno = tuple(
            chain.from_iterable(mol_anno.GetSubstructMatches(mcs_mol))
        )
        all_subs_ref = tuple(chain.from_iterable(mol_ref.GetSubstructMatches(mcs_mol)))
    else:
        all_subs_anno = ()
        all_subs_ref = ()

    ref_img = Draw.MolToImage(mol_ref, highlightAtoms=all_subs_ref, wedgeBonds=False)
    anno_img = Draw.MolToImage(mol_anno, highlightAtoms=all_subs_anno, wedgeBonds=False)
    return anno_img, ref_img


def show_ref_spectrums(spectrum_state, structure_obj, evt: gr.SelectData):
    try:
        line_num = evt.index[0]
        return get_reference_table(spectrum_state, structure_obj, line_num)
    except Exception as e:
        logging.error(f"Error in show_ref_spectrum: {e}")
        return None, None

def get_reference_table(spectrum_state, structure_obj, idx=0):
    try:                
        smi_anno = structure_obj["CanonicalSMILES"][idx]
        current_reference = spectrum_state.metadata["reference"]
        annotation = spectrum_state.metadata["annotation"]
        i = np.where(annotation["CanonicalSMILES"].values == smi_anno)[0][0]
        reference_table = []
        for s in current_reference:
            if "smiles" in s.metadata.keys():
                smiles = s.metadata["smiles"]
            else:
                smiles = ""
            if "compound_name" in s.metadata.keys():
                name = s.metadata["compound_name"]
            else:
                name = smiles
            if "adduct" in s.metadata.keys():
                adduct = s.metadata["adduct"]
            else:
                adduct = ""
            if "parent_mass" in s.metadata.keys():
                parent_mass = s.metadata["parent_mass"]
            else:
                parent_mass = ""
            if "database" in s.metadata.keys():
                ref_database = s.metadata["database"]
            else:
                ref_database = ""
            reference_table.append([name, adduct, smiles, parent_mass, ref_database])
        reference_table = pd.DataFrame(
            reference_table, columns=["name", "adduct", "smiles", "parent_mass", "database"]
        )  # 创建一个DataFrame对象，用于存储参考表格的数据

        return reference_table, smi_anno
    except Exception as e:
        logging.error(f"Error in get_reference_table: {e}")
        return pd.DataFrame(columns=["name", "adduct", "smiles", "parent_mass", "database"]), ""

def show_structure_select_all(spectrum_state, structure_obj, evt: gr.SelectData):
    try:
        # 先更新 Reference Spectrums 表格 & structure_state
        ref_spectrums, structure_state = show_ref_spectrums(spectrum_state, structure_obj, evt)
        # 再用刚更新的 structure_state 和 spectrum_state 直接画分子图
        anno_img, ref_img = show_default_mol(structure_state, spectrum_state)
        return ref_spectrums, structure_state, anno_img, ref_img
    except Exception as e:
        logging.error(f"Error in on_structure_select_all: {e}")
        return pd.DataFrame(columns=["name", "adduct", "smiles", "parent_mass", "database"]), "", None, None
def get_default_structure_select_all(spectrum_state, structure_obj):
    try:
        # 先更新 Reference Spectrums 表格 & structure_state
        ref_spectrums, structure_state = get_reference_table(spectrum_state, structure_obj)
        # 再用刚更新的 structure_state 和 spectrum_state 直接画分子图
        anno_img, ref_img = show_default_mol(structure_state, spectrum_state)
        return ref_spectrums, structure_state, anno_img, ref_img
    except Exception as e:
        logging.error(f"Error in on_structure_select_all: {e}")
        return pd.DataFrame(columns=["name", "adduct", "smiles", "parent_mass", "database"]), "", None, None
    
if __name__ == "__main__":
    print(dpi_config, width_config, length_config)
