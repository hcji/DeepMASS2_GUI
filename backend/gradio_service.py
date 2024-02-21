import hashlib
import os
import time
import zipfile
from zipfile import ZipFile

import gradio as gr
import numpy as np
import pandas as pd
from matchms.importing import load_from_mgf

from backend.utils.identify_unkown import id_spectrum_list
from backend.utils.plot_utils import get_formula_mass


def show_formula_from_state(res_state, idx):
    """
    从spectrum state中读取当前选中的spectrum
    Args:
        res_state: 当前选中的spectrum
        idx: 当前选中的idx

    Returns:

    """
    formula_list = []
    cur_spectrum = None
    try:
        formula_list = np.unique(
            res_state["Identified Spectrum"][idx].metadata["annotation"][
                "MolecularFormula"
            ]
        )
        cur_spectrum = res_state["Identified Spectrum"][idx]
    except:
        print("no identified spectrum found")

    mass = [get_formula_mass(f) for f in formula_list]

    if cur_spectrum and ("parent_mass" in cur_spectrum.metadata.keys()):
        diff = np.array(
            [abs(m - float(cur_spectrum.metadata["parent_mass"])) for m in mass]
        )
    else:
        diff = np.repeat(np.nan, len(mass))

    formula_df = pd.DataFrame(
        {"formula": formula_list, "mass": mass, "error (mDa)": 1000 * diff}
    )

    return cur_spectrum, formula_df


def show_structure(spectrum_state, evt: gr.SelectData):
    line_num = evt.index[0]
    formula_list = spectrum_state.metadata["annotation"]["MolecularFormula"]
    select_formula = formula_list[line_num]

    annotation = spectrum_state.metadata["annotation"]
    structural_table = annotation.loc[
        annotation["MolecularFormula"] == select_formula, :
    ]
    structural_table = structural_table.reset_index(drop=True)
    return select_formula, structural_table


def show_ref_spectrum(cur_spectrum, evt: gr.SelectData):
    line_num = evt.index[0]
    from backend.utils.plot_utils import plot_2_spectrum

    fig_loss = plot_2_spectrum(
        cur_spectrum, cur_spectrum.metadata["reference"][line_num], loss=True
    )
    fig = plot_2_spectrum(
        cur_spectrum, cur_spectrum.metadata["reference"][line_num], loss=False
    )
    return fig_loss, fig


def show_info(cur_spectrum, evt: gr.SelectData):
    line_num = evt.index[0]
    d = cur_spectrum.metadata["reference"][line_num].metadata
    df = pd.DataFrame.from_dict(d, orient="index", columns=["value"])
    df.reset_index(inplace=True)
    df.rename(columns={"index": "key"}, inplace=True)
    return df


def show_formula(res_state, evt: gr.SelectData):
    # print(evt)
    # print(evt.__dict__)
    print(f"You selected {evt.value} at {evt.index} from {evt.target}")

    # 从click事件中获取行号
    line_num = evt.index[0]

    return show_formula_from_state(res_state, line_num)


def load_files(file_list, request: gr.Request):
    # if os.path.getsize(file_list[0])>1024*1024:
    #     raise gr.Error('File size too large')
    if any(os.path.getsize(file) > 1024 * 1024 for file in file_list):
        raise gr.Error("File size too large")
    if len(file_list) > 500:
        raise gr.Error("Too many spectra, please upload a smaller number of files")

    try:
        spectrum_list = []
        for fileName in file_list:  # 遍历每一个文件名
            spectrum_list += [s for s in load_from_mgf(fileName)]
    except Exception as e:
        raise gr.Error("Please upload mgf file")

    titles = [
        s.metadata["compound_name"]
        if "compound_name" in list(s.metadata.keys())
        else f"spectrum {i}"
        for i, s in enumerate(spectrum_list)
    ]
    spectrums_df = pd.DataFrame({"title": titles, "spectrum": spectrum_list})
    # 用于返回nav的质谱名列表
    name_list = spectrums_df[["title"]]

    return spectrums_df, name_list


def deepms_click_fn(state_df, request: gr.Request, progress=gr.Progress()):
    """点击run deepms的按钮触发事件

    Args:
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
    try:
        from backend.entity.user import Work

        start_time = time.time()
        res = id_spectrum_list(state_df["spectrum"], progress)

    except Exception as e:
        raise gr.Error("Error processing data")

    # end_time = time.time()
    # work_duration = end_time - start_time
    # w = Work(
    #     contact_info=contract_info,
    #     spectrum_count=len(res),
    #     submit_time=start_time,
    #     end_time=end_time,
    #     work_duration=work_duration,
    # )
    # engine = create_engine("sqlite:///./User_Information.db", echo=True)
    # Session = sessionmaker(bind=engine)
    # session = Session()
    # session.add(w)
    # session.commit()

    state_df["Identified Spectrum"] = res
    (
        annotation,
        spectrum_state,
        formula_df,
        formula_state,
        structural_table,
        structure_state,
    ) = (None, None, None, None, None, None)

    try:
        if "annotation" in res[0].metadata.keys():
            annotation = res[0].metadata["annotation"]
        else:
            return state_df, spectrum_state, formula_state, structure_state, formula_df
        spectrum_state, formula_df = show_formula_from_state(state_df, 0)
        formula_state = annotation["MolecularFormula"][0]
        structural_table = annotation.loc[
            annotation["MolecularFormula"] == formula_state, :
        ]
        structure_state = structural_table["CanonicalSMILES"][0]
    except Exception as e:
        raise gr.Error(str(e))

    return state_df, spectrum_state, formula_state, structure_state, formula_df


def click_matchms_fn(state_df, progress=gr.Progress()):
    try:
        res = id_spectrum_list(state_df["spectrum"], progress, is_deepmass=False)

        state_df["Identified Spectrum"] = res

    except Exception as e:
        raise gr.Error("Error processing data")

    try:
        annotation = res[0].metadata["annotation"]
        spectrum_state, formula_df = show_formula_from_state(state_df, 0)

        formula_state = annotation["MolecularFormula"][0]
        structural_table = annotation.loc[
            annotation["MolecularFormula"] == formula_state, :
        ]
        structure_state = structural_table["CanonicalSMILES"][0]
        # TODO 拿到username，并向数据库中插入一条该用户的作业信息

    except:
        # raise gr.Error("Error processing data")
        spectrum_state, formula_state, structure_state, formula_df = (
            None,
            None,
            None,
            None,
        )
        pass

    return state_df, spectrum_state, formula_state, structure_state, formula_df


def save_identification_csv(res_state):
    try:
        file_list = []
        dir_path = "../temp"
        for s in res_state["Identified Spectrum"]:
            name = s.metadata["compound_name"]
            if "annotation" in s.metadata.keys():
                annotation = s.metadata["annotation"]
            else:
                annotation = pd.DataFrame(
                    columns=["Title", "MolecularFormula", "CanonicalSMILES", "InChIKey"]
                )
            path = os.path.join(dir_path, f"{name}.csv")
            csv = annotation.to_csv(path)
            file_list.append(path)
        md5_obj = hashlib.md5()
        md5_obj.update(str(file_list).encode("utf-8"))
        md5_name = md5_obj.hexdigest()
        zip_path = os.path.join(dir_path, f"{md5_name}.zip")
        with ZipFile(zip_path, "w") as zip_obj:
            for f in file_list:
                zip_obj.write(f, compress_type=zipfile.ZIP_DEFLATED)
        file_list.insert(0, zip_path)
        return gr.File(file_list, visible=True)
    except Exception as e:
        raise gr.Error("Error saving identification CSV")
