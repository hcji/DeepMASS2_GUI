#!/usr/bin/env bash
pwd
git remote add origin https://github.com/hcji/DeepMASS2_GUI
git reset --hard origin/lkr_dev
source activate deepmass2
python ./backend/gradio_app.py & python ./backend/register.py