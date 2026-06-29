#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
对 free_fitting_original.py 中拟合失败的两个星系（vagc_47886, vagc_13459）
进行专门拟合。

与原始版本的区别：
  - eps 固定为 0.1（原版 0.5）
  - initsma 从 sma_27 开始，逐步 /2, /3, /4, /5
  - 每个 initsma 都尝试 4 个 PA（0°, 45°, 90°, 135°）
  - 只处理指定的两个星系，不使用多进程
"""

from astropy.table import Table
from astropy.io import fits
import numpy as np
import os

from config import *
from my_tools import *


def run():
    """对 vagc_47886 和 vagc_13459 进行专门拟合"""

    # ---- 目标星系 ----
    targets = [
        ('vagc', 47886),
        ('vagc', 13459),
    ]

    # ---- 路径 ----
    dir_img  = IMG_DIR_2
    dir_mask = TOTAL_MASK_TOTAL
    dir_iso  = PROPS_ORIGINAL_ISOTAB

    # ---- 拟合参数 ----
    eps = 0.1
    pa_candidates = [0.0, np.pi/4, np.pi/2, 3*np.pi/4]  # 0°, 45°, 90°, 135°
    frac_list = [1, 2, 3, 4, 5]  # sma_27 / frac

    # ---- 读取总表 ----
    table = Table.read(TABLE_PATH)

    for label, index in targets:
        # 从表中筛选该星系
        mask = (table['survey'] == label) & (table['index'] == index)
        if not np.any(mask):
            print(f"{label}_{index} not found in table")
            continue

        row = table[mask][0]
        ra = row['ra']
        dec = row['dec']
        sma_27 = row['proc1_r27']

        # 文件路径
        suffix = f"{label}_{index}.fits"
        suffix_clr = f"a_band/{label}_{index}.fits"
        path_img  = os.path.join(dir_img,  suffix_clr)
        path_mask = os.path.join(dir_mask, suffix)
        path_iso  = os.path.join(dir_iso,  suffix)

        print(f"\n--- Processing {label}_{index} ---")
        print(f"  sma_27 = {sma_27:.2f} pix")

        # 读取图像 & 掩模
        img, hdr = fits.getdata(path_img, header=True)
        mask_data = fits.getdata(path_mask)
        masked_image = make_masked_image(img, mask_data)

        # 几何中心
        xc, yc = convert_wcs(ra, dec, hdr)
        minsma = 1.0
        maxsma = xc / 2

        # ---- 遍历 initsma × PA ----
        fitted = False
        for frac in frac_list:
            initsma = get_initsma(sma_27 / frac)
            print(f"  trying initsma = sma_27/{frac} = {initsma:.2f}")

            for pa in pa_candidates:
                pa_deg = pa * 180 / np.pi
                iso_table = try_fit(
                    masked_image, xc, yc, eps, pa,
                    initsma, minsma, maxsma,
                )
                if iso_table is not None:
                    iso_table.write(path_iso)
                    print(f"  ✓ SUCCESS: initsma={initsma:.2f}, pa={pa_deg:.0f}°")
                    fitted = True
                    break

            if fitted:
                break

        if not fitted:
            print(f"  ✗ {label}_{index} failed to fit")


if __name__ == "__main__":
    run()
