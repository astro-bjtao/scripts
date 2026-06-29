#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
根据mask_auto, mask_manual文件生成mask_total.
"""

from astropy.table import Table
from astropy.io import fits
import numpy as np
import os

# 路径配置
from config import *
# 通用工具
from my_tools import check_dir

def run_single(table, parms):
    """
    单线程处理函数：遍历星表，对每个星系合并掩模。
    """
    dir_mask_auto   = parms['mask_auto']
    dir_mask_manual = parms['mask_manual']
    dir_mask_total  = parms['mask_total']

    for i in range(len(table)):
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]

        suffix = f"{label}_{index}.fits"
        file_mask_auto   = os.path.join(dir_mask_auto, suffix)
        file_mask_manual = os.path.join(dir_mask_manual, suffix)
        out_file = os.path.join(dir_mask_total, suffix)

        # 读取内区和外区掩模
        mask_auto, header = fits.getdata(file_mask_auto, header=True)
        mask_manual = fits.getdata(file_mask_manual)

        # 合并
        mask_total = (mask_auto>0) | (mask_manual>0)

        # 保存
        fits.PrimaryHDU(data=mask_total.astype(np.uint8),
                        header=header).writeto(out_file)

def run_all():

    table = Table.read(TABLE_PATH)

    check_dir(TOTAL_MASK_TOTAL, clean=True)

    parms = {
        'mask_auto': TOTAL_MASK_AUTO,   # 来自 config.py
        'mask_manual': TOTAL_MASK_MANUAL,
        'mask_total': TOTAL_MASK_TOTAL
    }

    run_single(table, parms)

if __name__ == "__main__":
    run_all()