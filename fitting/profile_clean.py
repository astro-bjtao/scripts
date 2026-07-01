#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
在干净图像上重做 In_Ellipse 测光，输出新的等强度线表。

模式:
  python fitting/profile_clean.py        → image（干净科学图像）
  python fitting/profile_clean.py var    → var（干净方差图像）
"""

import os
from config import *
from my_tools import *


def process_chunk(tasks, parms):
    """处理 (label, index, clr, row_dict) 任务。"""
    dir_img = parms['img_dir']
    dir_iso = parms['iso_dir']
    dir_out = parms['out_dir']

    for label, index, clr, _ in tasks:
        suffix_clr = f"{clr}_band/{label}_{index}.fits"
        suffix     = f"{label}_{index}.fits"

        img = fits.getdata(os.path.join(dir_img, suffix_clr))
        iso_in = Table.read(os.path.join(dir_iso, suffix))
        iso = In_Ellipse(img, iso_in)
        iso_out = reshape_isotable(iso_to_table(iso))
        iso_out.write(os.path.join(dir_out, suffix_clr))


def run_image():
    for clr in ['g', 'r', 'i']:
        check_dir(os.path.join(PROPS_ORIGINAL_CLEAN_ISO, f"{clr}_band/"), clean=True)
    parms = {'img_dir': PROPS_ORIGINAL_CLEAN_IMAGE,
             'iso_dir': PROPS_ORIGINAL_ISOTAB,
             'out_dir': PROPS_ORIGINAL_CLEAN_ISO}
    run_multi_flat(process_chunk, TABLE_PATH, parms)


def run_var():
    for clr in ['g', 'r', 'i']:
        check_dir(os.path.join(PROPS_ORIGINAL_CLEAN_VAR_ISO, f"{clr}_band/"), clean=True)
    parms = {'img_dir': PROPS_ORIGINAL_CLEAN_VAR,
             'iso_dir': PROPS_ORIGINAL_ISOTAB,
             'out_dir': PROPS_ORIGINAL_CLEAN_VAR_ISO}
    run_multi_flat(process_chunk, TABLE_PATH, parms)


if __name__ == "__main__":
    import sys
    run_var() if sys.argv[-1] == 'var' else run_image()
