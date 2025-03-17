import base64
import os

import gradio as gr

from backend.service.gradio_service import (clear_files, deepms_click_fn,
                                            get_structure_data_frame,
                                            load_files,
                                            save_identification_csv,
                                            show_default_ref_spectrum,
                                            show_formula, show_ref_spectrum,
                                            show_structure)
from backend.utils.auth import auth_ps
from backend.utils.plot_utils import (get_reference_table, show_default_mol,
                                      show_mol, show_ref_spectrums)
from backend.utils.theme import Seafoam

# matplotlib.use('Agg')
os.environ["MPLCONFIGDIR"] = os.getcwd() + "/configs/"

seafoam = Seafoam()
import os

# path = "./frontend/icon/beian.png"
# print("绝对路径:", os.path.abspath(path))
# print("文件是否存在:", os.path.exists(path))

"""读取图片并转换为 base64 编码，返回 'data:image/png;base64,...' 形式的字符串"""
def get_beian_image_as_base64(img_path: str) -> str:
    
    with open(img_path, 'rb') as f:
        data = f.read()
    base64_str = base64.b64encode(data).decode('utf-8')
    return f"data:image/png;base64,{base64_str}"

with gr.Blocks(
    title="DeepMS 2", theme=seafoam, css="footer {visibility: hidden}"
) as demo:
    # 保存读取文件的结果
    # res_state = gr.Dataframe(visible=False)
    res_state = gr.State([])
    # 保存当前选择的spectrum
    spectrum_state = gr.State([])
    # 保存当前选择的structure
    structure_state = gr.State([])
    # 保存压缩文件目标名称
    target_zip_file_name_state = gr.State([])

    with gr.Row():
        file_obj = gr.File(file_count="multiple", type="filepath", height=100)
        download = gr.File(visible=False, interactive=False)
    with gr.Row(elem_classes=["first_row"]):
        run_save_btn = gr.Button("Save")
        run_deepms_btn = gr.Button(
            "Run DeepMS",
        )
        # run_matchms_btn = gr.Button("Run MatchMS")
    #
    with gr.Row():
        with gr.Column(scale=1):
            nav_obj = gr.DataFrame(
                headers=["name"],
                elem_classes=["scroll"],
                interactive=False,
                height=200,
                label="Navigator",
                wrap=True,
            )
        with gr.Column(scale=1):
            formula_obj = gr.DataFrame(
                headers=["Formula"],
                elem_classes=["scroll"],
                interactive=False,
                height=200,
                label="Formula Finder",
                wrap=True,
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
            wrap=True,
        )

    with gr.Row():
        ref_spectrums = gr.DataFrame(
            label="Reference Spectrums",
            headers=["name", "adduct", "smiles", "parent_mass", "database"],
            interactive=False,
            height=200,
            column_widths=["20%"],
            wrap=True,
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
            information_obj = gr.DataFrame(
                interactive=False,
                wrap=True,
            )

    # 上传文件自动更新
    file_obj.upload(
        load_files,
        inputs=file_obj,
        outputs=[res_state, nav_obj, target_zip_file_name_state],
    )
    # 删除文件清除状态
    file_obj.clear(
        clear_files,
        inputs=None,
        outputs=[
            res_state,
            spectrum_state,
            structure_state,
            nav_obj,
            formula_obj,
            structure_obj,
            ref_spectrums,
            spectrum_plot_fig,
            spectrum_loss_plot_fig,
            ann_structure_fig,
            ref_structure_fig,
            information_obj,
        ],
    )
    run_deepms_btn.click(
        fn=deepms_click_fn,
        inputs=[res_state],
        outputs=[res_state, spectrum_state, formula_obj, information_obj],
        concurrency_limit=4,
    )
    # run_matchms_btn.click(
    #     fn=click_matchms_fn,
    #     inputs=[res_state],
    #     outputs=[
    #         res_state,
    #         spectrum_state,
    #         formula_state,
    #         structure_state,
    #         formula_obj,
    #     ],
    # )

    # 选中Navigator的一行
    nav_obj.select(
        fn=show_formula,
        inputs=[res_state],
        outputs=[spectrum_state, formula_obj, information_obj],
    )
    # 选中Formula Finder的一行
    formula_obj.select(
        fn=show_structure,
        inputs=[
            spectrum_state,
        ],
        outputs=[structure_obj],
    )
    # Formula Finder的内容改变，，继续触发Structure Finder的内容改变
    formula_obj.change(
        fn=get_structure_data_frame,
        inputs=[
            spectrum_state,
        ],
        outputs=[structure_obj],
    )
    # 选中Structure Finder的一行
    structure_obj.select(
        fn=show_ref_spectrums,
        inputs=[spectrum_state, structure_obj],
        outputs=[ref_spectrums, structure_state],
    )
    # Structure Finder的内容改变，，继续触发Reference Spectrums的内容改变
    structure_obj.change(
        fn=get_reference_table,
        inputs=[spectrum_state, structure_obj],
        outputs=[ref_spectrums, structure_state],
    )

    # 选中Reference Spectrums的时候更新对比质谱图
    ref_spectrums.select(
        fn=show_ref_spectrum,
        inputs=[spectrum_state],
        outputs=[
            spectrum_loss_plot_fig,
            spectrum_plot_fig,
        ],
    )
    # 选中Reference Spectrums的时候更新对比质谱图
    ref_spectrums.change(
        fn=show_default_ref_spectrum,
        inputs=[spectrum_state],
        outputs=[
            spectrum_loss_plot_fig,
            spectrum_plot_fig,
        ],
    )
    # 选中Reference Spectrums的时候更新对比分子图
    ref_spectrums.select(
        fn=show_mol,
        inputs=[structure_state, spectrum_state],
        outputs=[
            ann_structure_fig,
            ref_structure_fig,
        ],
    )
    # 选中Reference Spectrums的时候更新对比分子图
    ref_spectrums.change(
        fn=show_default_mol,
        inputs=[structure_state, spectrum_state],
        outputs=[
            ann_structure_fig,
            ref_structure_fig,
        ],
    )
    # 点击保存按钮
    run_save_btn.click(
        fn=save_identification_csv,
        inputs=[res_state, target_zip_file_name_state],
        outputs=[download],
    )
    # 在所有组件之后添加备案信息组件
    beian_img_src = get_beian_image_as_base64("./frontend/icon/beian.png")
    # print("**************beian_img_src*************",beian_img_src[:100])
    gr.HTML(f"""
    <div style="width:360px; margin:36px auto 0 auto; 
                display: flex; 
                align-items: center; 
                justify-content: center; 
                gap: 16px; 
                font-size: 10px;">

        <a href="https://beian.mps.gov.cn/#/query/webSearch?code=43010402002088"
        target="_blank" rel="noreferrer"
        style="display: flex; align-items: center; text-decoration: none; color: #a5aaa3;">
        <img src="{beian_img_src}" alt="公安备案图标"
                style="width: 13px; height: 13px; margin-right: 4px;">
        <span>湘公网安备43010402002088号</span>
        </a>


        <a href="https://beian.miit.gov.cn/#/Integrated/index"
        target="_blank" rel="noreferrer"
        style="display: flex; align-items: center; text-decoration: none; color: #a5aaa3;">
        <span>湘ICP备2025102844号-1</span>
        </a>
    </div>
    """)


if __name__ == "__main__":
    print("Starting Webui!!!!")
    # demo.queue(concurrency_count=2)
    demo.launch(
        server_name="0.0.0.0",
        server_port=12341,
        show_api=False,
        auth=auth_ps,
        max_threads=2,
        auth_message="Please input your email and password",
    )
    print("Started Webui!!!!")
