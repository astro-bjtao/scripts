#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
合并 companion, segmap_all, star, unusable, 生成mask_total_auto。
"""

from astropy.table import Table
from astropy.io import fits
import numpy as np
from multiprocessing import Process
import time
import os

# 路径配置
from config import *
# 通用工具
from my_tools import check_dir


def run_single(table, parms):
    """
    单线程处理函数：遍历星表，对每个星系合并掩模。
    """
    dir_mask_companion  = parms['mask_companion']
    dir_mask_segmap_all = parms['mask_segmap_all']
    dir_mask_star       = parms['mask_star']
    dir_mask_unusable   = parms['mask_unusable']
    dir_mask_total_auto = parms['mask_total_auto']

    for i in range(len(table)):
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]

        suffix = f"{label}_{index}.fits"
        file_mask_companion  = os.path.join(dir_mask_companion, suffix)
        file_mask_segmap_all = os.path.join(dir_mask_segmap_all, suffix)
        file_mask_star       = os.path.join(dir_mask_star, suffix)
        file_mask_unusable   = os.path.join(dir_mask_unusable, suffix)
        out_file = os.path.join(dir_mask_total_auto, suffix)

        # 读取内区和外区掩模
        mask_companion, header = fits.getdata(file_mask_companion, header=True)
        mask_segmap_all = fits.getdata(file_mask_segmap_all)
        mask_star       = fits.getdata(file_mask_star)
        mask_unusable   = fits.getdata(file_mask_unusable)

        # 合并
        mask_total_auto = (mask_companion>0) | (mask_segmap_all>0) | (mask_star>0) | (mask_unusable>0)

        # 保存
        fits.PrimaryHDU(data=mask_total_auto.astype(np.uint8),
                        header=header).writeto(out_file)


def run_multi(path_table, parms, ncpu=120):
    """
    多进程调度：拆表 → 分配进程 → 合并。
    """
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
    """
    主入口：设置路径，准备输出目录，启动多进程处理。
    """
    check_dir(TOTAL_MASK_AUTO, clean=True)

    parms = {
        'mask_companion': COMPANION_MASK_DIR,   # 来自 config.py
        'mask_segmap_all': MASK_SEGMAP_ALL_AUTO,
        'mask_star': STAR_MASK_DIR,
        'mask_unusable': UNUSABLE_MASK,
        'mask_total_auto': TOTAL_MASK_AUTO
    }

    run_multi(TABLE_PATH, parms, ncpu=120)


if __name__ == "__main__":
    run_all()