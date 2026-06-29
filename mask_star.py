#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成亮星掩模 (star mask)，g/r/i 波段并集，排除中心目标区域。
阈值：star_model > 0.0178 视为亮星区域。
输出：二值掩模 (uint8)，保存至 Process2/cutout_mask/mask_star/mask/。
"""

from astropy.table import Table
from astropy.io import fits
import numpy as np
from multiprocessing import Process
import time
import os

# 架构导入
from config import *
from my_tools import check_dir

def run_single(table, parms):
    dir_mdl     = parms['star_mdl']
    dir_target  = parms['mask_target']
    dir_mask    = parms['star_mask']

    for i in range(len(table)):
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]

        # 读取 mask_target（中心区域）
        path_target = os.path.join(dir_target, f"{label}_{index}.fits")
        try:
            target, hdr = fits.getdata(path_target, header=True)
        except FileNotFoundError:
            print(f"[WARNING] Missing target mask for {label}_{index}, skipping.")
            continue

        mask_union = None
        for clr in ['g', 'r', 'i']:
            path_mdl = os.path.join(dir_mdl, f"{clr}_band/{label}_{index}.fits")
            try:
                img, _ = fits.getdata(path_mdl, header=True)
                if mask_union is None:
                    mask_union = np.zeros_like(img, dtype=bool)
                mask_union = mask_union | (img > 0.02)
            except FileNotFoundError:
                print(f"[WARNING] Missing {clr}-band star model for {label}_{index}")

        if mask_union is None:
            print(f"[ERROR] No star model found for {label}_{index}, skipping.")
            continue

        # 排除中心目标区域
        mask_star = mask_union & (target == 0)

        out_file = os.path.join(dir_mask, f"{label}_{index}.fits")
        fits.PrimaryHDU(data=mask_star.astype(np.uint8),
                        header=hdr).writeto(out_file)

def run_multi(path_table, parms, ncpu=120):
    tab = Table.read(path_table)

    n_total = len(tab)
    indices = np.linspace(0, n_total, ncpu + 1, dtype=int)
    process_list = []
    for j in range(ncpu):
        start = indices[j]
        end = indices[j + 1]
        if start == end:
            continue
        sub_tab = tab[start:end]
        process_list.append(
            Process(target=run_single, args=(sub_tab, parms))
        )

    for p in process_list:
        p.start()
        time.sleep(0.01)

    for p in process_list:
        p.join()

def run_all():
    check_dir(STAR_MASK_DIR, clean=True)

    parms = {
        'star_mdl':    STAR_MODEL_DIR,
        'mask_target': MASK_30_TARGET,
        'star_mask':   STAR_MASK_DIR
    }

    run_multi(TABLE_PATH, parms, ncpu=120)

if __name__ == "__main__":
    run_all()