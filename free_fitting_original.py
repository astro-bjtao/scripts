#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
对减过背景的星系图像进行测光。
"""

from astropy.table import Table, Column
from astropy.io import fits
import numpy as np
import shutil, os

# 路径配置
from config import *
# 通用工具
from my_tools import *

### fitting and show result

def run_single(table, parms):

    dir_img  = parms['image_dir']
    dir_mask = parms['mask_dir']
    dir_iso  = parms['isotab_dir']

    # pa 候选值（弧度）：0°, 45°, 90°, 135°
    pa_candidates = [0.0, np.pi/4, np.pi/2, 3*np.pi/4]

    for i in range(len(table)):
        # information
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]
        ra    = table['ra'].data[i]
        dec   = table['dec'].data[i]
        sma_27   = table['proc1_r27'].data[i]
        # paths
        suffix = f"{label}_{index}.fits"
        suffix_clr = f"a_band/{label}_{index}.fits"
        path_img     = os.path.join(dir_img,   suffix_clr)
        path_mask    = os.path.join(dir_mask,  suffix)
        path_iso     = os.path.join(dir_iso,   suffix)
        # image
        img, hdr = fits.getdata(path_img, header=True)
        mask = fits.getdata(path_mask)
        masked_image = make_masked_image(img, mask)
        # geometry
        xc, yc = convert_wcs(ra, dec, hdr)
        minsma = 1.0
        maxsma = xc/2
        eps = 0.5

        # 几何序列，保证不同星系的 sma 对齐
        list_init = (300/0.168) * 1.1**(-1*np.arange(60))

        def get_initsma(target):
            """从 list_init 中找最接近 target 的值"""
            delta = np.abs(list_init - target)
            flag = delta < (np.min(delta) + 0.01)
            return float(list_init[flag][0])

        def try_fit(initsma, pa):
            """尝试拟合，成功返回 iso_table，失败返回 None"""
            try:
                iso = Ellipse_free(masked_image,
                                    xc, yc, eps, pa,
                                    initsma, minsma, maxsma,
                                    step=0.1, fix_center=True, fix_pa=False)
                iso_tab = iso_to_table(iso)
                iso_table = reshape_isotable(iso_tab)
                if len(iso_table) > 0:
                    return iso_table
            except Exception:
                pass
            return None

        # 阶段 1：逐步减小 initsma（pa=0）
        fitted = False
        for frac in [2, 3, 4, 5]:
            initsma = get_initsma(sma_27 / frac)
            iso_table = try_fit(initsma, pa_candidates[0])
            if iso_table is not None:
                iso_table.write(path_iso)
                fitted = True
                break

        # 阶段 2：固定 initsma=sma_27/5，尝试不同 pa
        if not fitted:
            initsma = get_initsma(sma_27 / 5)
            for pa in pa_candidates[1:]:  # pa=0 已在阶段 1 尝试过
                iso_table = try_fit(initsma, pa)
                if iso_table is not None:
                    iso_table.write(path_iso)
                    fitted = True
                    break

        if not fitted:
            print(f"{label}_{index} failed to fit")

def run_all():
    """
    主入口：设置路径，准备输出目录，启动多线程处理。
    """
    check_dir(PROPS_ORIGINAL_ISOTAB, clean=True)

    parms = {
        'image_dir':   IMG_DIR_2,
        'mask_dir':    TOTAL_MASK_TOTAL,
        'isotab_dir':  PROPS_ORIGINAL_ISOTAB      
    }

    run_multi(run_single, TABLE_PATH, parms, ncpu=120)
if __name__ == "__main__":

    run_all()