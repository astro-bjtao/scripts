#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
自动生成外区 (outer) 掩模，保留 segmap_15 的标签值。

逻辑：
    1. 以 segmap_15 的标签图作为初始化掩模。
    2. 遍历 segmap_05 中的每一个非零标签。
    3. 对每个标签内的像素，进行连通域分析。
    4. 若某连通域与 segmap_15 有重叠，且不与 MASK_15_TARGET 有任何交集，
       则找出重叠区域内出现次数最多的 segmap_15 标签值，将该连通域所有
       像素赋值为该标签。
    5. 所采用的按标签孤立再分连通域的做法，杜绝了不同标签在空间上被
       错误合并的可能。
    6. 若连通域触及 MASK_15_TARGET（中心目标区域），则直接丢弃，防止
       掩模污染星系自身。
"""

from astropy.table import Table
from astropy.io import fits
import numpy as np
from scipy.ndimage import label as scilabel
import os
import shutil

# 路径配置
from config import *
# 通用工具
from my_tools import check_dir, run_multi


def run_single(table, parms):
    """
    单线程处理函数：遍历星表，对每个星系生成外区标记掩模。
    """
    dir_mask_05 = parms['mask_05']
    dir_mask_15 = parms['mask_15']
    dir_mask_target = parms['mask_target']   # 新增：mask_15 的目标区域
    dir_mask_out = parms['mask_out']

    for i in range(len(table)):
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]

        suffix = f"{label}_{index}.fits"
        file_mask_05 = os.path.join(dir_mask_05, suffix)
        file_mask_15 = os.path.join(dir_mask_15, suffix)
        file_mask_target = os.path.join(dir_mask_target, suffix)  # 新增
        out_file = os.path.join(dir_mask_out, suffix)

        # 读取 segmap_15、segmap_05 和目标区域
        mask_15, header = fits.getdata(file_mask_15, header=True)
        mask_05 = fits.getdata(file_mask_05)
        mask_target = fits.getdata(file_mask_target)  # 新增

        # 初始化输出掩模：直接继承 segmap_15 的标签值
        mask = mask_15.copy()

        # 1. 获取 mask_05 中所有唯一的非零标签
        unique_labels = np.unique(mask_05)
        unique_labels = unique_labels[unique_labels > 0]

        # 2. 遍历每一个标签
        for lab in unique_labels:
            binary = (mask_05 == lab)

            if not np.any(binary):
                continue

            # 3. 在当前标签内进行连通域分析
            area_labels, N = scilabel(binary.astype(np.int32))

            # 4. 遍历每一个连通域
            for k in range(1, N + 1):
                flag_05 = (area_labels == k)

                # 5. 关键检查：若连通域涉及目标区域，则丢弃
                if np.any(flag_05 & (mask_target > 0)):
                    continue

                # 6. 检查与 mask_15 是否有重叠
                overlap_pixels = mask_15[flag_05]
                overlap_values = overlap_pixels[overlap_pixels > 0]

                if len(overlap_values) > 0:
                    # 取重叠区域内出现次数最多的标签（众数）
                    assigned_label = np.bincount(overlap_values).argmax()
                    # 将整个连通域赋值为该标签
                    mask[flag_05] = assigned_label

        # 保存最终的外区标签掩模
        fits.PrimaryHDU(data=mask, header=header).writeto(out_file)


def run_all():
    """
    主入口：设置路径，准备输出目录，启动多进程处理，复制 mask_target。
    """
    # 输出目录
    check_dir(MASK_OUTER_AUTO, clean=True)

    parms = {
        'mask_05': MASK_05,
        'mask_15': MASK_15,
        'mask_target': MASK_15_TARGET,   # 新增参数
        'mask_out': MASK_OUTER_AUTO
    }

    run_multi(run_single, TABLE_PATH, parms, ncpu=120)

    # 复制 mask_target（从 segmap_15 的 mask_target 到 mask_outer 下）
    check_dir(MASK_OUTER_TARGET, clean=True)
    for item in os.listdir(MASK_15_TARGET):
        src = os.path.join(MASK_15_TARGET, item)
        dst = os.path.join(MASK_OUTER_TARGET, item)
        if os.path.isfile(src):
            shutil.copy2(src, dst)


if __name__ == "__main__":
    run_all()