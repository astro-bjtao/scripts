#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
用椭圆模型替换被掩模像素，生成干净图像。

支持两种模式：
  - image : 科学图像（IMG_DIR_2）+ bmodel → clean_image/
  - var   : 方差图像（VAR_DIR）+ bmodel_var → clean_var/
"""

import numpy as np
import os
from astropy.io import fits

from config import *
from my_tools import *


def run_single(table, parms):
    dir_img    = parms['img_dir']
    dir_mask   = parms['mask_dir']
    dir_model  = parms['model_dir']
    dir_out    = parms['out_dir']
    is_var     = parms.get('is_var', False)

    for i in range(len(table)):
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]

        for clr in ['g', 'r', 'i']:
            suffix_clr = f"{clr}_band/{label}_{index}.fits"
            suffix     = f"{label}_{index}.fits"

            path_img   = os.path.join(dir_img,   suffix_clr)
            path_mask  = os.path.join(dir_mask,  suffix)
            path_model = os.path.join(dir_model, suffix_clr)
            path_out   = os.path.join(dir_out,   suffix_clr)

            img, hdr = fits.getdata(path_img, header=True)
            if is_var:
                img[np.isnan(img) | np.isinf(img)] = 0.0

            mask  = fits.getdata(path_mask)
            model = fits.getdata(path_model)

            # 用模型替换被掩模像素
            img[mask > 0] = model[mask > 0]
            img = img.astype(np.float32)

            fits.PrimaryHDU(data=img, header=hdr).writeto(path_out, overwrite=True)


def run_image():
    """科学图像 clean → PROPS_ORIGINAL_CLEAN_IMAGE"""
    for clr in ['g', 'r', 'i']:
        check_dir(os.path.join(PROPS_ORIGINAL_CLEAN_IMAGE, f"{clr}_band/"), clean=True)

    parms = {
        'img_dir':   IMG_DIR_2,
        'mask_dir':  TOTAL_MASK_TOTAL,
        'model_dir': PROPS_ORIGINAL_BMODEL,
        'out_dir':   PROPS_ORIGINAL_CLEAN_IMAGE,
    }
    run_multi(run_single, TABLE_PATH, parms, ncpu=120)


def run_var():
    """误差图像 clean → PROPS_ORIGINAL_CLEAN_VAR"""
    for clr in ['g', 'r', 'i']:
        check_dir(os.path.join(PROPS_ORIGINAL_CLEAN_VAR, f"{clr}_band/"), clean=True)

    parms = {
        'img_dir':   VAR_DIR,
        'mask_dir':  TOTAL_MASK_TOTAL,
        'model_dir': PROPS_ORIGINAL_BMODEL_VAR,
        'out_dir':   PROPS_ORIGINAL_CLEAN_VAR,
        'is_var':    True,
    }
    run_multi(run_single, TABLE_PATH, parms, ncpu=120)


if __name__ == "__main__":
    import sys
    mode = sys.argv[1] if len(sys.argv) > 1 else 'image'
    if mode == 'var':
        run_var()
    else:
        run_image()
