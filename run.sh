#!/usr/bin/env bash
pwd
source activate deepmass2
python ./backend/gradio_app.py & python ./backend/register.py