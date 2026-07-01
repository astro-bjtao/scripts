#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
用椭圆模型替换被掩模像素，生成干净图像。

模式:
  python fitting/clean_pixels.py        → image（科学图像）
  python fitting/clean_pixels.py var    → var（方差图像）

并行: 星系×波段展开，均分到 120 进程。
"""

import numpy as np
import os
from astropy.io import fits

from config import *
from my_tools import *


def process_chunk(tasks, parms):
    """处理一批 (label, index, clr, row_dict) 任务。"""
    dir_img   = parms['img_dir']
    dir_mask  = parms['mask_dir']
    dir_model = parms['model_dir']
    dir_out   = parms['out_dir']
    is_var    = parms.get('is_var', False)

    for label, index, clr, _ in tasks:
        suffix_clr = f"{clr}_band/{label}_{index}.fits"
        suffix     = f"{label}_{index}.fits"

        img, hdr = fits.getdata(os.path.join(dir_img,   suffix_clr), header=True)
        if is_var:
            img[np.isnan(img) | np.isinf(img)] = 0.0

        mask  = fits.getdata(os.path.join(dir_mask,  suffix))
        model = fits.getdata(os.path.join(dir_model, suffix_clr))

        img[mask > 0] = model[mask > 0]
        img = img.astype(np.float32)

        fits.PrimaryHDU(data=img, header=hdr).writeto(
            os.path.join(dir_out, suffix_clr), overwrite=True)


def run_image():
    for clr in ['g', 'r', 'i']:
        check_dir(os.path.join(PROPS_ORIGINAL_CLEAN_IMAGE, f"{clr}_band/"), clean=True)
    parms = {'img_dir': IMG_DIR_2, 'mask_dir': TOTAL_MASK_TOTAL,
             'model_dir': PROPS_ORIGINAL_BMODEL, 'out_dir': PROPS_ORIGINAL_CLEAN_IMAGE}
    run_multi_flat(process_chunk, TABLE_PATH, parms)


def run_var():
    for clr in ['g', 'r', 'i']:
        check_dir(os.path.join(PROPS_ORIGINAL_CLEAN_VAR, f"{clr}_band/"), clean=True)
    parms = {'img_dir': VAR_DIR, 'mask_dir': TOTAL_MASK_TOTAL,
             'model_dir': PROPS_ORIGINAL_BMODEL_VAR, 'out_dir': PROPS_ORIGINAL_CLEAN_VAR,
             'is_var': True}
    run_multi_flat(process_chunk, TABLE_PATH, parms)


if __name__ == "__main__":
    import sys
    mode = sys.argv[1] if len(sys.argv) > 1 else 'image'
    run_var() if mode == 'var' else run_image()
