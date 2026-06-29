#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成伴星系 (companion) 掩模。

逻辑：
    对于每个目标星系，读取其 SExtractor 星表和分割图，
    找出面积大于目标面积一半 或 MAG_AUTO < 20 的明亮伴星系，
    将这些伴星系周围 2.5 * R90 的椭圆区域标记为掩模。
    输出：标签掩模 (uint8)，保存至 Process2/cutout_mask/mask_companion/。
"""

from astropy.table import Table
from astropy.io import fits, ascii
import numpy as np
import os

# 架构导入
from config import *
from my_tools import check_dir, get_ellipse, run_multi

def run_single(table, parms):
    """
    对星表中的每行生成伴星系掩模。
    """
    dir_cat    = parms['cat_dir']
    dir_seg    = parms['seg_dir']
    dir_target = parms['target_dir']
    dir_out    = parms['mask_dir']

    for i in range(len(table)):
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]
        # R27 参数（用于计算目标面积）
        r27_sma  = table['proc1_r27'].data[i]  # 半长轴
        r27_eps  = 1.0 - table['q_disk'].data[i] if 'q_disk' in table.colnames else 0.5  # 椭率近似
        # 若星表中没有 EPS，暂时用 0.5 代替；你可以根据实际情况调整
        
        # 文件路径
        path_cat    = os.path.join(dir_cat,    f"{label}_{index}.cat")
        path_seg    = os.path.join(dir_seg,    f"{label}_{index}.fits")
        path_target = os.path.join(dir_target, f"{label}_{index}.fits")
        path_out    = os.path.join(dir_out,    f"{label}_{index}.fits")

        # 检查必要文件
        if not (os.path.exists(path_cat) and os.path.exists(path_seg) and os.path.exists(path_target)):
            print(f"[WARNING] Missing input files for {label}_{index}, skipping.")
            continue

        # 读取数据
        cat    = ascii.read(path_cat)
        seg, hdr = fits.getdata(path_seg, header=True)
        target = fits.getdata(path_target)

        # ---------- 核心逻辑：扩展伴星系 ----------
        # 目标星系面积
        b2a_target = 1.0 - r27_eps
        area_target = np.pi * b2a_target * r27_sma * r27_sma
        area_companion = area_target / 2.0

        # 找出在目标区域内的 seg 标签
        in_target = target > 0
        in_target_labels = np.unique(seg[in_target])
        in_target_labels = in_target_labels[in_target_labels > 0]

        # 可扩展的伴星系条件
        crit_1 = cat['ISOAREA_IMAGE'].data > area_companion
        crit_2 = cat['MAG_AUTO'].data < 20
        available_indices = np.where(crit_1 | crit_2)[0]

        # 生成扩展分割图
        ext_seg = np.zeros(seg.shape, dtype=np.int32)
        for j in available_indices:
            if (j + 1) in in_target_labels:  # 跳过目标本身
                continue

            cenx = cat['X_IMAGE'].data[j] - 1
            ceny = cat['Y_IMAGE'].data[j] - 1
            r90  = cat['FLUX_RADIUS'].data[j]
            sma  = cat['A_IMAGE'].data[j]
            smb  = cat['B_IMAGE'].data[j]
            b2a  = smb / sma + 0.1
            if b2a > 1.0:
                b2a = 1.0
            pa   = cat['THETA_IMAGE'].data[j] - 90

            ellipse = get_ellipse(seg.shape, cenx, ceny, b2a, pa)
            in_ellipse = ellipse < (r90 * 2.5)
            ext_seg[in_ellipse] = j + 1

        # 移除目标区域
        ext_seg[target > 0] = 0

        # 保存
        fits.PrimaryHDU(data=ext_seg,
                        header=hdr).writeto(path_out)

def run_all():
    check_dir(COMPANION_MASK_DIR, clean=True)

    parms = {
        'cat_dir':    COMPANION_CAT_DIR,
        'seg_dir':    COMPANION_SEG_DIR,
        'target_dir': COMPANION_TARGET_DIR,
        'mask_dir':   COMPANION_MASK_DIR
    }

    run_multi(run_single, TABLE_PATH, parms, ncpu=120)

if __name__ == "__main__":
    run_all()