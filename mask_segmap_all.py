#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
合并 inner 和 outer 掩模，生成最终的 mask_all。

逻辑：
    1. 分别读取 mask_inner 和 mask_outer。
    2. 对 mask_inner 膨胀 2 像素，mask_outer 膨胀 3 像素（填补掩模间隙）。
    3. 取并集，生成 mask_all。
    4. 由于 mask_outer 保留标签值而 mask_inner 是二值图，
       当两者重叠时优先保留 mask_outer 的标签（众数），
       仅在 mask_outer 为空但 mask_inner 有值时使用统一标记值 999。
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
from my_tools import check_dir, extend_mask


def run_single(table, parms):
    """
    单线程处理函数：遍历星表，对每个星系合并掩模。
    """
    dir_mask_inner = parms['mask_inner']
    dir_mask_outer = parms['mask_outer']
    dir_mask_all = parms['mask_all']

    for i in range(len(table)):
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]

        suffix = f"{label}_{index}.fits"
        file_mask_inner = os.path.join(dir_mask_inner, suffix)
        file_mask_outer = os.path.join(dir_mask_outer, suffix)
        out_file = os.path.join(dir_mask_all, suffix)

        # 读取内区和外区掩模
        mask_inner, header = fits.getdata(file_mask_inner, header=True)
        mask_outer = fits.getdata(file_mask_outer)

        # 膨胀操作（填补掩模边缘间隙）
        mask_inner = extend_mask(mask_inner, ext_r=2)
        mask_outer = extend_mask(mask_outer, ext_r=3)

        # 合并策略：优先保留 mask_outer 的标签值
        mask_all = np.zeros_like(mask_outer, dtype=np.uint8)

        # 条件1：mask_outer > 0 → 直接继承标签
        mask_all[mask_outer > 0] = mask_outer[mask_outer > 0]

        # 条件2：mask_outer == 0 但 mask_inner > 0 → 用 2 标记
        flag_only_inner = (mask_outer == 0) & (mask_inner > 0)
        mask_all[flag_only_inner] = 2

        # 保存
        fits.PrimaryHDU(data=mask_all.astype(np.uint8),
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
    check_dir(MASK_SEGMAP_ALL_AUTO, clean=True)

    parms = {
        'mask_inner': MASK_INNER_AUTO,   # 来自 config.py
        'mask_outer': MASK_OUTER_AUTO,
        'mask_all': MASK_SEGMAP_ALL_AUTO
    }

    run_multi(TABLE_PATH, parms, ncpu=120)


if __name__ == "__main__":
    run_all()