import logging
import os
import uuid
import zipfile
from zipfile import ZipFile

import gradio as gr
import numpy as np
import pandas as pd

from backend.load_config import GLOBAL_CONFIG
from backend.service.job_service import JobService
from backend.utils.common_utils import get_title_from_spectrum
from backend.utils.identify_unkown import id_spectrum_list
from backend.utils.plot_utils import get_formula_mass
from backend.utils.plot_utils import plot_2_spectrum
from backend.utils.spectrum_process import load_spectrum_file

MAX_SPECTRUM_NUM = GLOBAL_CONFIG["identification"]["max_spectrum"]
MAX_FILES = GLOBAL_CONFIG["identification"]["max_files"]
MAX_FILE_SIZE = eval(GLOBAL_CONFIG["identification"]["max_file_size"])


def validate_spectrum_annotation(spectrum):
    pass


def show_formula_from_state(res_state, idx):
    """
    从spectrum state中读取当前选中的spectrum
    Args:
        res_state: 当前选中的spectrum
        idx: 当前选中的idx

    Returns:当前选中质谱
            当前选中质谱的候选分子式表格Dataframe
            当前选中质谱的元数据

    """
    formula_list = []
    # 判断是否已经被鉴定
    if "Identified Spectrum" not in res_state.keys():
        raise gr.Error("Missing Identified Spectrum")
    # 选中质谱
    cur_spectrum = res_state["Identified Spectrum"][idx]
    # 判断当前选中质谱是否拥有注释
    if "annotation" not in cur_spectrum.metadata.keys():
        raise gr.Error("Missing Identified Spectrum annotation")
    annotation = cur_spectrum.metadata["annotation"]
    # 判断返回结果是否含有分子式
    if "MolecularFormula" not in annotation.keys():
        raise gr.Error("Missing Annotation MolecularFormula")
    # 获取分子式列表
    formula_list = np.unique(annotation["MolecularFormula"])
    # 获取分子式的相对分子质量
    mass = [get_formula_mass(f) for f in formula_list]
    # 获取相对分子质量与母离子质量差
    if cur_spectrum and ("parent_mass" in cur_spectrum.metadata.keys()):
        parent_mass = float(cur_spectrum.metadata["parent_mass"])
        diff = np.array([abs(m - parent_mass) for m in mass])
    else:
        diff = np.repeat(np.nan, len(mass))
        gr.Warning("Spectrum Missing Parent Mass")
    # 分子式表格
    formula_df = pd.DataFrame(
        {"formula": formula_list, "mass": mass, "error (mDa)": 1000 * diff}
    )
    formula_df = formula_df.round(3)

    return cur_spectrum, formula_df, show_info(cur_spectrum)


def show_structure(spectrum_state, evt: gr.SelectData):
    """
    Show the structure table
    展示候选结构表
    于点击分子式时触发
    Args:
        spectrum_state: 当前选中的质谱
        evt: 点击事件

    Returns:选择的分子式
            候选结构表

    """
    # 获取点击行号
    line_num = evt.index[0]

    structural_table = get_structure_data_frame(spectrum_state, line_num)
    return structural_table


def get_structure_data_frame(spectrum_state, idx=0):
    formula_list = spectrum_state.metadata["annotation"]["MolecularFormula"]
    select_formula = formula_list[idx]

    annotation = spectrum_state.metadata["annotation"]
    structural_table = annotation.loc[
        annotation["MolecularFormula"] == select_formula, :
    ]
    structural_table = structural_table.reset_index(drop=True)
    structural_table = structural_table.round(3)
    return structural_table


def show_ref_spectrum(cur_spectrum, evt: gr.SelectData):
    line_num = evt.index[0]
    return show_default_ref_spectrum(cur_spectrum, line_num)


def show_default_ref_spectrum(cur_spectrum, idx=0):
    fig_loss = plot_2_spectrum(
        cur_spectrum, cur_spectrum.metadata["reference"][idx], loss=True
    )
    fig = plot_2_spectrum(
        cur_spectrum, cur_spectrum.metadata["reference"][idx], loss=False
    )
    return fig_loss, fig


def show_info(cur_spectrum):
    logging.info("点击结构行")
    if "reference" not in cur_spectrum.metadata.keys():
        gr.Error("No reference")
    try:
        d = cur_spectrum.metadata.copy()
        del d["reference"]
        del d["annotation"]
    except:
        logging.info("没有选中的质谱")
        return pd.DataFrame()
    # 输出质谱元数据
    df = pd.DataFrame.from_dict(d, orient="index", columns=["value"])
    df.reset_index(inplace=True)
    df.rename(columns={"index": "key"}, inplace=True)
    logging.info("获取质谱元数据成功")
    return df


def show_formula(res_state, evt: gr.SelectData):
    """

    Args:
        res_state:鉴定结果
        evt:点击事件

    Returns: 当前选中的质谱，当前质谱候选分子式的DataFrame

    """
    print(
        f"Spectrum List Click----selected {evt.value} at {evt.index} from {evt.target}"
    )

    # 从click事件中获取行号
    line_num = evt.index[0]

    return show_formula_from_state(res_state, line_num)


def load_files(file_list, request: gr.Request):
    """
    从文件列表中读取质谱，并输出到窗口
    Args:
        file_list: 文件名列表，是保存在服务器临时路径中
        request: Gradio自带Reqyest对象，可用于获取信息

    Returns:
        所有读取的质谱
        所有质谱的显示名称
    """
    is_super_user = request.username == "admin@deepmass.cn"
    # 判断文件大小总和，不接受大于一定阈值的请求
    sum_file_size = 0
    for file in file_list:
        sum_file_size += os.path.getsize(file)
    if sum_file_size > MAX_FILE_SIZE and (not is_super_user):
        raise gr.Error("File size too large")
    # 目标保存压缩文件名称
    target_zip_file_name = "result.zip"
    if len(file_list) == 1:
        target_zip_file_name = f"{os.path.basename(file_list[0])}.zip"
    # 读取每个质谱文件
    spectrum_list = []
    for file_name in file_list:  # 遍历每一个文件名
        try:
            loaded_spectra_list = load_spectrum_file(file_name)
        except Exception as e:
            raise gr.Error("Please upload standard file")
        spectrum_list.extend(loaded_spectra_list)
        if len(spectrum_list) > MAX_SPECTRUM_NUM and (not is_super_user):
            raise gr.Error(
                f"Only a maximum of {MAX_SPECTRUM_NUM} spectra are allowed to be uploaded"
            )

    # 获取所有的质谱名称，若无，则使用编号代替
    titles = [
        s.metadata["compound_name"]
        if "compound_name" in list(s.metadata.keys())
        else f"spectrum {i}"
        for i, s in enumerate(spectrum_list)
    ]
    spectrums_df = pd.DataFrame({"title": titles, "spectrum": spectrum_list})
    # 用于返回nav的质谱名列表
    name_list = spectrums_df[["title"]]

    return spectrums_df, name_list, target_zip_file_name


def deepms_click_fn(state_df, request: gr.Request, progress=gr.Progress()):
    """点击run deepms的按钮触发事件

    Args:
        progress:
        request:
        state_df (_type_): _description_
        输入为一个dataframe,列名为title,spectrum

    Returns:
        _type_: _description_
        更新下列状态
            res_state,增加 identified spectrum字段,内为注释过的spectrum对象
            spectrum_state,设置选中的spectrum
            formula_state,,设置选中的formula
            structure_state,,设置选中的structure

    """
    # print(request.__dict__)
    contract_info = request.username
    logging.info(f"用户{contract_info}点击了run deepmass")
    try:
        job_id = JobService().start_job(contract_info)
        gr.Info("Start Identifing")
        res = id_spectrum_list(state_df["spectrum"], progress)
        JobService().end_job(job_id)
    except Exception as e:
        logging.error(e)
        raise gr.Error("Error processing data")

    state_df["Identified Spectrum"] = res
    gr.Info("Identified Successed")
    # spectrum_state = res[0] if len(res) > 0 else None
    cur_spectrum, formula_df, info_df = get_default_formula_dataframe(state_df)
    return state_df, cur_spectrum, formula_df, info_df


def get_default_formula_dataframe(state_df):
    return show_formula_from_state(state_df, 0)


# def click_matchms_fn(state_df, progress=gr.Progress()):
#     try:
#         res = id_spectrum_list(state_df["spectrum"], progress, is_deepmass=False)
#
#         state_df["Identified Spectrum"] = res
#
#     except Exception as e:
#         raise gr.Error("Error processing data")
#
#     try:
#         annotation = res[0].metadata["annotation"]
#         spectrum_state, formula_df = show_formula_from_state(state_df, 0)
#
#         formula_state = annotation["MolecularFormula"][0]
#         structural_table = annotation.loc[
#             annotation["MolecularFormula"] == formula_state, :
#         ]
#         structure_state = structural_table["CanonicalSMILES"][0]
#
#     except:
#         # raise gr.Error("Error processing data")
#         spectrum_state, formula_state, structure_state, formula_df = (
#             None,
#             None,
#             None,
#             None,
#         )
#         pass
#
#     return state_df, spectrum_state, formula_state, structure_state, formula_df


def save_identification_csv(res_state, target_zip_file_name_state):
    gr.Info('Saving identification CSV"')
    file_list = []
    dir_path = "./backend/temp/"
    dir_path = os.path.join(dir_path, uuid.uuid4().hex)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    # 判断是否有鉴定结果
    if res_state is None or "Identified Spectrum" not in res_state.columns:
        gr.Error("Missing Identified results, Run Identify First")
    for idx, s in enumerate(res_state["Identified Spectrum"]):
        # name = s.metadata["compound_name"]
        name = get_title_from_spectrum(spectrum=s, idx=idx)
        # 判断是否有鉴定结果，若无，则给空表
        if "annotation" in s.metadata.keys():
            annotation = s.metadata["annotation"]
        else:
            annotation = pd.DataFrame(
                columns=["Title", "MolecularFormula", "CanonicalSMILES", "InChIKey"]
            )
        # 生成文件名
        path = os.path.join(dir_path, f"{name}.csv")
        # 转csv文件
        csv = annotation.to_csv(path)
        file_list.append(path)
    # 计算压缩包文件名和路径
    zip_path = os.path.join(dir_path, target_zip_file_name_state)
    # 生成压缩文件
    with ZipFile(zip_path, "w") as zip_obj:
        for f in file_list:
            zip_obj.write(
                f, compress_type=zipfile.ZIP_DEFLATED, arcname=os.path.basename(f)
            )
    file_list.insert(0, zip_path)
    # 返回文件列表，第一个是压缩包，其他的是单个谱的鉴定结果
    return gr.File(file_list, visible=True)


def clear_files():
    """
    删除文件时，清空所有显示的状态
    Returns:

    """
    return None, None, None, None, None, None, None, None, None, None, None, None
