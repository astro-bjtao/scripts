#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
根据等强度线构建椭圆模型图像。

支持两种模式：
  - image : 对科学图像（IMG_DIR_2）建模 → bmodel/
  - var   : 对误差图像（VAR_DIR，sigma→variance）建模 → bmodel_var/
"""

import numpy as np
import os
from astropy.io import fits
from photutils.isophote import build_ellipse_model

from config import *
from my_tools import *


def run_single(table, parms):
    dir_img    = parms['img_dir']
    dir_mask   = parms['mask_dir']
    dir_iso    = parms['iso_dir']
    dir_out    = parms['out_dir']
    is_var     = parms.get('is_var', False)

    for i in range(len(table)):
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]

        for clr in ['g', 'r', 'i']:
            suffix_clr = f"{clr}_band/{label}_{index}.fits"
            suffix     = f"{label}_{index}.fits"

            path_img  = os.path.join(dir_img,  suffix_clr)
            path_mask = os.path.join(dir_mask, suffix)
            path_iso  = os.path.join(dir_iso,  suffix)
            path_out  = os.path.join(dir_out,  suffix_clr)

            img, hdr = fits.getdata(path_img, header=True)

            if is_var:
                # sigma → variance
                img[np.isnan(img) | np.isinf(img)] = 0.0
                img = img ** 2

            mask = fits.getdata(path_mask)
            masked = make_masked_image(img, mask)

            iso_table = Table.read(path_iso)
            iso = In_Ellipse(masked, iso_table)
            model = build_ellipse_model(img.shape, iso)
            model = model.astype(np.float32)

            fits.PrimaryHDU(data=model, header=hdr).writeto(path_out, overwrite=True)


def run_image():
    """科学图像椭圆建模 → PROPS_ORIGINAL_BMODEL"""
    for clr in ['g', 'r', 'i']:
        check_dir(os.path.join(PROPS_ORIGINAL_BMODEL, f"{clr}_band/"), clean=True)

    parms = {
        'img_dir':  IMG_DIR_2,
        'mask_dir': TOTAL_MASK_TOTAL,
        'iso_dir':  PROPS_ORIGINAL_ISOTAB,
        'out_dir':  PROPS_ORIGINAL_BMODEL,
    }
    run_multi(run_single, TABLE_PATH, parms, ncpu=120)


def run_var():
    """误差图像椭圆建模 → PROPS_ORIGINAL_BMODEL_VAR"""
    for clr in ['g', 'r', 'i']:
        check_dir(os.path.join(PROPS_ORIGINAL_BMODEL_VAR, f"{clr}_band/"), clean=True)

    parms = {
        'img_dir':  VAR_DIR,
        'mask_dir': TOTAL_MASK_TOTAL,
        'iso_dir':  PROPS_ORIGINAL_ISOTAB,
        'out_dir':  PROPS_ORIGINAL_BMODEL_VAR,
        'is_var':   True,
    }
    run_multi(run_single, TABLE_PATH, parms, ncpu=120)


if __name__ == "__main__":
    import sys
    mode = sys.argv[1] if len(sys.argv) > 1 else 'image'
    if mode == 'var':
        run_var()
    else:
        run_image()
