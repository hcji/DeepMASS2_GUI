import gradio as gr
from gradio.themes import Base


class Seafoam(Base):
    def __init__(self):
        super().__init__()
        self.name = "Seafoam"
        super().set(
            background_fill_primary="*primary_200",
            shadow_drop="0 1px 4px 0 rgb(240 248 255 / 0.5)",
            shadow_drop_lg="0 2px 5px 0 rgb(0 197 205 / 0.7)",
            block_background_fill="*primary_100",
            block_label_padding="*spacing_sm *spacing_md",
            block_label_background_fill="*primary_100",
            block_label_radius="*radius_md",
            block_label_text_size="*text_md",
            block_label_text_weight="600",
            block_label_text_color="*primary_400",
            button_shadow="*shadow_drop_lg",
            button_primary_background_fill="*primary_300",
            button_primary_text_color="black",
            button_secondary_background_fill="PowderBlue",
            button_secondary_text_color="*black",
            block_border_width="1px",
            panel_border_width="1px",
        )


def evaluate(username, issue, category, improvement, contact):
    if not issue:  # 确保问题字段是必填的
        return "请填写存在的问题。"
    feedback = f"收到来自{'匿名' if not username else username}的反馈：\n问题：{issue}\n问题类别：{category}\n改进建议：{improvement}\n联系方式：{'未提供' if not contact else contact}"
    return feedback


seafoam_theme = Seafoam()
username_input = gr.Textbox(label="用户名（选填）")
issue_input = gr.Textbox(lines=2, placeholder="请输入存在的问题（必填）", label="存在的问题")
category_input = gr.CheckboxGroup(choices=["技术问题", "服务问题", "其他"], label="问题类别")
improvement_input = gr.Textbox(label="改进建议（选填）")
contact_input = gr.Textbox(label="联系方式（选填）")

iface = gr.Interface(
    fn=evaluate,
    inputs=[
        username_input,
        issue_input,
        category_input,
        improvement_input,
        contact_input,
    ],
    outputs="text",
    title="用户反馈系统",
    # description="请填写您对系统的评价。",
    theme=seafoam_theme,
    css="footer {display: none !important;}",
)

iface.launch()
