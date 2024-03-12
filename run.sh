#!/usr/bin/env bash
pwd
git remote add origin https://github.com/hcji/DeepMASS2_GUI
git reset --hard origin/lkr_dev
source activate deepmass2
lsof -i:12341 | grep 'TCP' | awk '{print $2}' | xargs kill -9
lsof -i:8000 | grep 'TCP' | awk '{print $2}' | xargs kill -9
python ./backend/gradio_app.py & python ./backend/register.py