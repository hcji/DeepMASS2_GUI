import os

os.environ["MPLCONFIGDIR"] = os.getcwd() + "/configs/"
from itertools import chain
from hnswlib import Index
from gensim.models import Word2Vec
from core.identification import identify_unknown, match_spectrum

# matplotlib.use('Agg')
from matplotlib.figure import Figure
import gradio as gr
from matchms.importing import load_from_mgf
import numpy as np
import pandas as pd
import os
import pickle
from matchms import filtering as msfilters
from matchms import Spectrum
from rdkit import Chem
from rdkit.Chem import Draw, rdFMCS
from molmass import Formula
from zipfile import ZipFile
import hashlib
import zipfile

import time
from gradio.themes.base import Base


class Seafoam(Base):
    def __init__(
        self,
    ):
        super().__init__()
        self.name = "Seafoam"
        super().set(
            # Colors
            background_fill_primary="*primary_50",
            slider_color="*primary_500",
            slider_color_dark="*primary_600",
            # Shadows
            shadow_drop="0 1px 4px 0 rgb(240 248 255 / 0.5)",
            shadow_drop_lg="0 2px 5px 0 rgb(0 197 205 / 0.7)",
            # Block Labels
            block_background_fill="	AliceBlue",
            block_label_padding="*spacing_sm *spacing_md",
            block_label_background_fill="*primary_100",
            block_label_background_fill_dark="*primary_600",
            block_label_radius="*radius_md",
            block_label_text_size="*text_md",
            block_label_text_weight="600",
            block_label_text_color="*primary_400",
            block_label_text_color_dark="white",
            block_title_radius="*block_label_radius",
            block_title_padding="*block_label_padding",
            block_title_background_fill="*block_label_background_fill",
            block_title_text_weight="600",
            block_title_text_color="*primary_500",
            block_title_text_color_dark="white",
            block_label_margin="*spacing_md",
            # Buttons
            shadow_spread="6px",
            button_shadow="*shadow_drop_lg",
            button_shadow_hover="*shadow_drop_lg",
            checkbox_label_shadow="*shadow_drop_lg",
            button_shadow_active="*shadow_inset",
            button_primary_background_fill="*primary_500",
            button_primary_background_fill_hover="*primary_100",
            button_primary_background_fill_hover_dark="*primary_500",
            button_primary_text_color="white",
            button_secondary_background_fill="		LightCyan",
            button_secondary_background_fill_hover="*neutral_100",
            button_secondary_background_fill_hover_dark="*primary_100",
            button_secondary_text_color="*neutral_800",
            # Borders
            block_border_width="0px",
            panel_border_width="1px",
        )


seafoam = Seafoam()

default_index_positive = "data/references_index_positive_spec2vec.bin"
default_index_negative = "data/references_index_negative_spec2vec.bin"
default_reference_positive = "data/references_spectrums_positive.pickle"
default_reference_negative = "data/references_spectrums_negative.pickle"
print("Start Loading database")
default_database = pd.read_csv("data/database.csv")
print("Start Loading Word2Vec")
deepmass_positive = Word2Vec.load("model/Ms2Vec_allGNPSpositive.hdf5")
deepmass_negative = Word2Vec.load("model/Ms2Vec_allGNPSnegative.hdf5")
print("Start Loading negative reference")

with open(default_reference_negative, "rb") as file:
    reference_negative = pickle.load(file)
print("Start Loading positive reference")
with open(default_reference_positive, "rb") as file:
    reference_positive = pickle.load(file)
print("Start Loading hnsw index")
index_negative = Index(space="l2", dim=300)
index_negative.load_index(default_index_negative)

index_positive = Index(space="l2", dim=300)
index_positive.load_index(default_index_positive)

precursors_positive = np.array([s.get("precursor_mz") for s in reference_positive])
precursors_negative = np.array([s.get("precursor_mz") for s in reference_negative])
print("Finish!!!")
print("-" * 100)


def get_formula_mass(formula):
    f = Formula(formula)
    return f.isotope.mass


def pre_process_spectrum(spectrum: Spectrum):
    spectrum = msfilters.add_parent_mass(spectrum)
    msfilters.add_precursor_mz(spectrum)
    m


def identify_pos(spectrum):
    return identify_unknown(
        spectrum,
        index_positive,
        deepmass_positive,
        reference_positive,
        default_database,
    )


def identify_neg(spectrum):
    return identify_unknown(
        spectrum,
        index_negative,
        deepmass_negative,
        reference_negative,
        default_database,
    )


def match_pos(spectrum):
    return match_spectrum(spectrum, precursors_positive, reference_positive)


def match_neg(spectrum):
    return match_spectrum(spectrum, precursors_negative, reference_negative)


def plot_2_spectrum(spectrum, reference, loss=False):
    mz, abunds = spectrum.peaks.mz, spectrum.peaks.intensities
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
            mz, abunds = spectrum.losses.mz, spectrum.losses.intensities
            mz1, abunds1 = reference.losses.mz, reference.losses.intensities
        except:
            print("Cannot Plot Losses")
            return
            abunds /= np.max(abunds)
    abunds1 /= np.max(abunds1)

    fig = Figure(figsize=(2, 1), dpi=300)
    fig.subplots_adjust(top=0.95, bottom=0.3, left=0.18, right=0.95)

    axes = fig.add_subplot(111)
    axes.tick_params(width=0.8, labelsize=3)
    axes.spines["bottom"].set_linewidth(0.5)
    axes.spines["left"].set_linewidth(0.5)
    axes.spines["right"].set_linewidth(0.5)
    axes.spines["top"].set_linewidth(0.5)
    axes.tick_params(width=0.8, labelsize=3)
    axes.vlines(mz, ymin=0, ymax=abunds, color="r", lw=0.5)
    axes.vlines(mz1, ymin=0, ymax=-abunds1, color="b", lw=0.5)
    axes.axhline(y=0, color="black", lw=0.5)
    axes.set_xlabel("m/z", fontsize=3.5)
    axes.set_ylabel("abundance", fontsize=3.5)
    return fig


def show_ref_spectrum(cur_spectrum, evt: gr.SelectData):
    line_num = evt.index[0]
    fig_loss = plot_2_spectrum(
        cur_spectrum, cur_spectrum.metadata["reference"][line_num], loss=True
    )
    fig = plot_2_spectrum(
        cur_spectrum, cur_spectrum.metadata["reference"][line_num], loss=False
    )
    return fig_loss, fig


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
        allsubs_anno = tuple(chain.from_iterable(mol_anno.GetSubstructMatches(mcs_mol)))
        allsubs_ref = tuple(chain.from_iterable(mol_ref.GetSubstructMatches(mcs_mol)))
    else:
        allsubs_anno = ()
        allsubs_ref = ()

    ref_img = Draw.MolToImage(mol_ref, highlightAtoms=allsubs_ref, wedgeBonds=False)
    anno_img = Draw.MolToImage(mol_anno, highlightAtoms=allsubs_anno, wedgeBonds=False)
    return anno_img, ref_img


def show_mol(structure_state, cur_spectrum, evt: gr.SelectData):
    line_num = evt.index[0]
    ref_smi = cur_spectrum.metadata["reference"][line_num].metadata["smiles"]
    anno_img, ref_img = plot_2_mol(structure_state, ref_smi)
    return anno_img, ref_img


def show_info(cur_spectrum, evt: gr.SelectData):
    line_num = evt.index[0]
    d = cur_spectrum.metadata["reference"][line_num].metadata
    df = pd.DataFrame.from_dict(d, orient="index", columns=["value"])
    df.reset_index(inplace=True)
    df.rename(columns={"index": "key"}, inplace=True)
    return df


def show_ref_spectrums(spectrum_state, structure_obj, evt: gr.SelectData):
    line_num = evt.index[0]
    smi_anno = structure_obj["CanonicalSMILES"][line_num]
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


def show_formula_from_state(res_state, idx):
    formula_list = np.unique(
        res_state["Identified Spectrum"][idx].metadata["annotation"]["MolecularFormula"]
    )
    cur_spectrum = res_state["Identified Spectrum"][idx]
    mass = [get_formula_mass(f) for f in formula_list]

    if "parent_mass" in cur_spectrum.metadata.keys():
        diff = np.array(
            [abs(m - float(cur_spectrum.metadata["parent_mass"])) for m in mass]
        )
    else:
        diff = np.repeat(np.nan, len(mass))

    formula_df = pd.DataFrame(
        {"formula": formula_list, "mass": mass, "error (mDa)": 1000 * diff}
    )

    return cur_spectrum, formula_df


def show_formula(res_state, evt: gr.SelectData):
    # print(evt)
    # print(evt.__dict__)
    print(f"You selected {evt.value} at {evt.index} from {evt.target}")

    # 从click事件中获取行号
    line_num = evt.index[0]

    return show_formula_from_state(res_state, line_num)


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


def id_spectrum_list(spectrum_list, progress=None, is_deepmass=True):
    res = []
    if is_deepmass:
        for s in progress.tqdm(spectrum_list):
            sn = None
            try:
                if "ionmode" in s.metadata.keys():
                    if s.metadata["ionmode"] == "negative":
                        sn = identify_neg(s)
                    else:
                        sn = identify_pos(s)
                else:
                    sn = identify_pos(s)
            except:
                pass
            res.append(sn)
    else:
        for s in progress.tqdm(spectrum_list):
            sn = None
            try:
                if "ionmode" in s.metadata.keys():
                    if s.metadata["ionmode"] == "negative":
                        sn = match_neg(s)
                    else:
                        sn = match_pos(s)
                else:
                    sn = match_pos(s)
            except:
                pass
            res.append(sn)
    return res


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
        from user import Work

        start_time = time.time()
        res = id_spectrum_list(state_df["spectrum"], progress)

    except Exception as e:
        raise gr.Error("Error processing data")

    end_time = time.time()
    work_duration = end_time - start_time
    w = Work(
        contact_info=contract_info,
        spectrum_count=len(res),
        submit_time=start_time,
        end_time=end_time,
        work_duration=work_duration,
    )
    engine = create_engine("sqlite:///./User_Information.db", echo=True)
    Session = sessionmaker(bind=engine)
    session = Session()
    session.add(w)
    session.commit()

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

        annotation = res[0].metadata["annotation"]
        spectrum_state, formula_df = show_formula_from_state(state_df, 0)

        formula_state = annotation["MolecularFormula"][0]
        structural_table = annotation.loc[
            annotation["MolecularFormula"] == formula_state, :
        ]
        structure_state = structural_table["CanonicalSMILES"][0]
        # TODO 拿到username，并向数据库中插入一条该用户的作业信息

        return state_df, spectrum_state, formula_state, structure_state, formula_df
    except Exception as e:
        raise gr.Error("Error processing data")


def save_identification_csv(res_state):
    try:
        file_list = []
        dir_path = "./temp"
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


with gr.Blocks(
    title="DeepMS 2", theme=seafoam, css="footer {visibility: hidden}"
) as demo:
    # 保存读取文件的结果
    # res_state = gr.Dataframe(visible=False)
    res_state = gr.State([])
    # 保存当前选择的spectrum
    spectrum_state = gr.State([])
    # 保存当前选择的formula
    formula_state = gr.State([])
    # 保存当前选择的structure
    structure_state = gr.State([])

    with gr.Row(elem_classes=["first_row"]):
        file_obj = gr.File(file_count="multiple", type="filepath", height=100)
        download = gr.File(visible=False, interactive=False)
    with gr.Row(elem_classes=["first_row"]):
        run_save_btn = gr.Button("Save")
        run_deepms_btn = gr.Button(
            "Run DeepMS",
        )
        run_matchms_btn = gr.Button("Run MatchMS")
    #
    with gr.Row(elem_classes=["secend_row"]):
        with gr.Column(scale=1):
            nav_obj = gr.DataFrame(
                headers=["name"],
                elem_classes=["scroll"],
                interactive=False,
                height=200,
                label="Navigator",
            )
        with gr.Column(scale=1):
            formula_obj = gr.DataFrame(
                headers=["Formula"],
                elem_classes=["scroll"],
                interactive=False,
                height=200,
                label="Formula Finder",
            )
    with gr.Row():
        structure_obj = gr.DataFrame(
            headers=[
                "Title",
                "MolecularFormula",
                "CanonicalSMILES",
                "InChIKey",
                "DeepMass Score",
            ],
            interactive=False,
            elem_classes=["scroll"],
            height=200,
            label="Structure Finder",
        )

    with gr.Row():
        ref_spectrums = gr.DataFrame(
            label="Reference Spectrums",
            headers=["name", "adduct", "smiles", "parent_mass", "database"],
            interactive=False,
            height=200,
            column_widths="20%",
        )
    with gr.Row():
        with gr.Tab(label="Spectrum"):
            with gr.Row():
                spectrum_plot_fig = gr.Plot(label="Spectrum")
                spectrum_loss_plot_fig = gr.Plot(label="Loss")
        with gr.Tab(label="Structure"):
            with gr.Row():
                ann_structure_fig = gr.Image(
                    label="Annotated Structure", height=200, width=200
                )
                ref_structure_fig = gr.Image(
                    label="Reference Structure", height=200, width=200
                )
        with gr.Tab(label="Information"):
            information_obj = gr.DataFrame(interactive=False)

    # 上传文件自动更新
    file_obj.change(
        load_files,
        inputs=file_obj,
        outputs=[
            res_state,
            nav_obj,
        ],
    )

    nav_obj.select(
        fn=show_formula, inputs=[res_state], outputs=[spectrum_state, formula_obj]
    )
    formula_obj.select(
        fn=show_structure,
        inputs=[
            spectrum_state,
        ],
        outputs=[formula_state, structure_obj],
    )
    structure_obj.select(
        fn=show_ref_spectrums,
        inputs=[spectrum_state, structure_obj],
        outputs=[ref_spectrums, structure_state],
    )

    run_deepms_btn.click(
        fn=deepms_click_fn,
        inputs=[res_state],
        outputs=[
            res_state,
            spectrum_state,
            formula_state,
            structure_state,
            formula_obj,
        ],
        concurrency_limit=4,
    )
    run_matchms_btn.click(
        fn=click_matchms_fn,
        inputs=[res_state],
        outputs=[
            res_state,
            spectrum_state,
            formula_state,
            structure_state,
            formula_obj,
        ],
    )

    ref_spectrums.select(
        fn=show_ref_spectrum,
        inputs=[spectrum_state],
        outputs=[
            spectrum_loss_plot_fig,
            spectrum_plot_fig,
        ],
    )

    ref_spectrums.select(
        fn=show_mol,
        inputs=[structure_state, spectrum_state],
        outputs=[
            ann_structure_fig,
            ref_structure_fig,
        ],
    )
    ref_spectrums.select(
        fn=show_info, inputs=[spectrum_state], outputs=[information_obj]
    )
    run_save_btn.click(
        fn=save_identification_csv, inputs=[res_state], outputs=[download]
    )

if __name__ == "__main__":
    print("Starting Webui!!!!")
    # demo.queue(concurrency_count=2)
    from sqlalchemy import select
    from sqlalchemy.orm import sessionmaker
    from sqlalchemy import create_engine

    from user import User
    from user import Login

    # TODO
    def auth_ps(contact_info, passwd):
        try:
            engine = create_engine("sqlite:///./User_Information.db", echo=True)
            Session = sessionmaker(bind=engine)
            session = Session()
            # TODO 在数据库中查询username与passwd是否匹配，username设置为unique
            passwd = hashlib.new("md5", passwd.encode("utf-8")).hexdigest()
            results = session.execute(
                select(User)
                .where(User.contact_info.in_([contact_info]))
                .where(User.passwd.in_([passwd]))
            )
            res = results.scalars().all()
            login_time = time.time()
            log = Login(contact_info=contact_info, login_time=login_time)
            session.add(log)
            session.commit()

            if len(res) > 0:
                return True
            else:
                return False
        except Exception as e:
            raise gr.Error("Error in authentication process")

    demo.launch(
        server_name="0.0.0.0",
        server_port=12341,
        show_api=False,
        auth=auth_ps,
        max_threads=2,
        auth_message="Please input your email and password",
    )
    print("Started Webui!!!!")
