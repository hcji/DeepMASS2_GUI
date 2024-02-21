from gradio.themes import Base


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
