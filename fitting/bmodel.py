#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
根据拟合的等强度线构建椭圆模型图像。

读取 isotab → 重建 Ellipse 几何 → build_ellipse_model → 输出模型 FITS。
"""

from photutils.isophote import build_ellipse_model

from astropy.io import fits
import numpy as np
import os

from config import *
from my_tools import *


def get_model(shape, isophotes):
    """从 IsophoteList 构建椭圆模型图像。"""
    return build_ellipse_model(shape, isophotes)


def run_single(table, parms):
    dir_img    = parms['img_dir']
    dir_mask   = parms['mask_dir']
    dir_iso    = parms['iso_dir']
    dir_bmodel = parms['bmodel_dir']

    for i in range(len(table)):
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]

        for clr in ['g', 'r', 'i']:
            suffix_clr = f"{clr}_band/{label}_{index}.fits"
            suffix     = f"{label}_{index}.fits"

            path_img    = os.path.join(dir_img,    suffix_clr)
            path_mask   = os.path.join(dir_mask,   suffix)
            path_iso    = os.path.join(dir_iso,    suffix)
            path_bmodel = os.path.join(dir_bmodel, suffix_clr)

            # 读取图像和掩模
            img, hdr = fits.getdata(path_img, header=True)
            mask = fits.getdata(path_mask)
            masked_image = make_masked_image(img, mask)

            # 读取等强度线，重建几何，构建模型
            iso_table = Table.read(path_iso)
            iso = In_Ellipse(masked_image, iso_table)
            model = get_model(img.shape, iso)
            model = model.astype(np.float32)

            fits.PrimaryHDU(data=model, header=hdr).writeto(
                path_bmodel, overwrite=True)


def run_all():
    """主入口。"""
    for clr in ['g', 'r', 'i']:
        check_dir(os.path.join(PROPS_ORIGINAL_BMODEL, f"{clr}_band/"), clean=True)

    parms = {
        'img_dir':    IMG_DIR_2,
        'mask_dir':   TOTAL_MASK_TOTAL,
        'iso_dir':    PROPS_ORIGINAL_ISOTAB,
        'bmodel_dir': PROPS_ORIGINAL_BMODEL,
    }

    run_multi(run_single, TABLE_PATH, parms, ncpu=120)


if __name__ == "__main__":
    run_all()
