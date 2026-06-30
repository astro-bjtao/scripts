#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
对拟合失败的星系，手动指定初始条件重新拟合。

初始条件根据目视检查设置：
  - eps=0.7（目视估计椭圆率）
  - pa 根据图像方向手动设定
  - initsma=r27/2
"""

from astropy.table import Table
from astropy.io import fits
import numpy as np
import os

from config import *
from my_tools import *


def run():
    """对手动指定初始条件的星系进行拟合"""

    # ---- 目标星系 (label, index, pa_deg_Xplus, eps) ----
    # PA 用 moments_estimate 估计（X+ 约定），initsma = r27/3
    targets = [
        # ('gama', 2573,  170, 0.7),   # 已跑通
        ('vagc', 47858, 82,  0.60),
        # ('vagc', 95552, 113, 0.55),  # 已跑通
        # ('vagc', 99948, 30,  0.8),   # 已跑通
    ]

    # initsma 分母：1 → sma_27，2 → r27/2，3 → r27/3
    initsma_frac = 1

    # ---- 路径 ----
    dir_img  = IMG_DIR_2
    dir_mask = TOTAL_MASK_TOTAL
    dir_iso  = PROPS_ORIGINAL_ISOTAB

    # ---- 读取总表 ----
    table = Table.read(TABLE_PATH)

    for label, index, pa_deg, eps in targets:
        # 从表中筛选该星系
        mask = (table['survey'] == label) & (table['index'] == index)
        if not np.any(mask):
            print(f"{label}_{index} not found in table")
            continue

        row = table[mask][0]
        ra    = row['ra']
        dec   = row['dec']
        r27   = row['proc1_r27']

        # PA 转换为弧度
        pa_rad = pa_deg * np.pi / 180.0

        # initsma 从几何序列中取最接近 r27/N 的值
        initsma = get_initsma(r27 / initsma_frac)

        # 文件路径
        suffix     = f"{label}_{index}.fits"
        suffix_clr = f"a_band/{label}_{index}.fits"
        path_img   = os.path.join(dir_img,  suffix_clr)
        path_mask  = os.path.join(dir_mask, suffix)
        path_iso   = os.path.join(dir_iso,  suffix)

        print(f"\n--- Processing {label}_{index} ---")
        print(f"  r27 = {r27:.2f}  pix")
        print(f"  eps = {eps},  pa = {pa_deg}°,  initsma = {initsma:.2f}")

        # 读取图像 & 掩模
        img, hdr = fits.getdata(path_img, header=True)
        mask_data = fits.getdata(path_mask)
        masked_image = make_masked_image(img, mask_data)

        # 几何中心
        xc, yc = convert_wcs(ra, dec, hdr)
        minsma = 1.0
        maxsma = xc / 2

        # 拟合
        iso_table = try_fit(
            masked_image, xc, yc, eps, pa_rad,
            initsma, minsma, maxsma,
        )

        if iso_table is not None:
            iso_table.write(path_iso, overwrite=True)
            print(f"  ✓ SUCCESS")
        else:
            print(f"  ✗ failed to fit")


if __name__ == "__main__":
    run()
