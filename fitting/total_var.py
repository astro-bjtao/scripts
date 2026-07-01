#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
计算总方差：模型方差 + 背景方差 + 低频噪声方差。

背景方差 = (bkg_valerr × area)²
低频方差 = (err_low / √(area/area_quarter))² × area²

依赖表列: bkg_lfqerr_{clr}, bkg_valerr_{clr}, bkg_rmin, bkg_rmax
"""

import numpy as np
import os
from config import *
from my_tools import *


def compute_total_var(sma, eps, tflux_e_var, bkg_valerr, err_low, bkg_r):
    """计算总方差数组。"""
    area = np.pi * sma * sma * (1.0 - eps)

    # 背景减除引入的方差
    back_var = (bkg_valerr * area) ** 2

    # 低频噪声：用背景环区域面积归一化
    area_quarter = (np.pi * (bkg_r * 1.1) ** 2 - np.pi * bkg_r ** 2) / 4.0
    num_region = np.maximum(area / area_quarter, 1.0)
    err_low_norm = err_low / np.sqrt(num_region)
    low_var = (err_low_norm * area) ** 2

    return tflux_e_var + back_var + low_var


def process_chunk(tasks, parms):
    for label, index, clr, row in tasks:
        suffix_clr = f"{clr}_band/{label}_{index}.fits"
        path_in  = os.path.join(parms['iso_dir'], suffix_clr)
        path_out = os.path.join(parms['out_dir'], suffix_clr)

        err_low    = row[f'bkg_lfqerr_{clr}']
        bkg_valerr = row[f'bkg_valerr_{clr}']
        bkg_r = np.sqrt(row['bkg_rmin'] * row['bkg_rmax'])

        iso = Table.read(path_in)
        iso['tflux_e_var_tot'] = compute_total_var(
            iso['sma'].data, iso['ellipticity'].data,
            iso['tflux_e'].data, bkg_valerr, err_low, bkg_r)
        iso.write(path_out)


def run():
    for clr in ['g', 'r', 'i']:
        check_dir(os.path.join(PROPS_ORIGINAL_VAR_TOTAL, f"{clr}_band/"), clean=True)
    parms = {'iso_dir': PROPS_ORIGINAL_CLEAN_VAR_ISO,
             'out_dir': PROPS_ORIGINAL_VAR_TOTAL}
    run_multi_flat(process_chunk, TABLE_PATH, parms)


if __name__ == "__main__":
    run()
