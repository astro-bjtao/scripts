#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
根据.reg文件生成mask.
"""

from astropy.table import Table
from astropy.io import fits
import numpy as np
import os

# 路径配置
from config import *
# 通用工具
from my_tools import check_dir, reg_to_mask

def run_single(table, parms):
    """
    单线程处理函数：遍历星表，对每个星系合并掩模。
    """
    dir_mask_total  = parms['mask_total']
    dir_mask_manual = parms['mask_manual']
    dir_reg_ds9     = parms['reg_ds9']

    for i in range(len(table)):
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]

        suffix = f"{label}_{index}.fits"
        file_mask_total  = os.path.join(dir_mask_total, suffix)
        out_file = os.path.join(dir_mask_manual, suffix)

        suffix = f"{label}_{index}.reg"
        file_reg_ds9 = os.path.join(dir_reg_ds9, suffix)

        try:
            reg_to_mask(file_reg_ds9, file_mask_total, out_file)
        except:
            print(f"{label}_{index} have a bad ds9 mask")

def run_all():

    table = Table.read(TABLE_PATH)

    check_dir(TOTAL_MASK_MANUAL, clean=True)

    parms = {
        'mask_total': TOTAL_MASK_AUTO,   # 来自 config.py
        'mask_manual': TOTAL_MASK_MANUAL,
        'reg_ds9': TOTAL_REG_DS9
    }

    run_single(table, parms)

if __name__ == "__main__":
    run_all()