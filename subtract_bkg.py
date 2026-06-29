#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
减去背景值。
"""

from astropy.table import Table
from astropy.io import fits
import numpy as np

import shutil
import os

# 路径配置
from config import *
# 通用工具
from my_tools import check_dir, MeasureBkg, run_multi

def run_single(table, parms):

    # parms
    dir_img  = parms['image_dir']
    dir_iso  = parms['isotab_dir']
    dir_outimg  = parms['image_out_dir']

    for i in range(len(table)):
        # information
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]
        bkg_rmin = table['bkg_rmin'].data[i]
        # path
        for clr in ['g', 'r', 'i']:
            # paths
            suffix_clr = f"{clr}_band/{label}_{index}.fits"
            path_img    = os.path.join(dir_img,    suffix_clr)
            path_iso    = os.path.join(dir_iso,    suffix_clr)
            path_outimg = os.path.join(dir_outimg, suffix_clr)
            # get local bkg
            isotab = Table.read(path_iso)
            rms, err_high, err_low, err_tot, local_bkg = MeasureBkg(isotab, bkg_rmin)
            # sub local bkg
            img, hdr = fits.getdata(path_img, header=True)
            img = img - local_bkg
            # save
            fits.PrimaryHDU(data=img, header=hdr).writeto(path_outimg)

def run_all():
    """
    主入口：设置路径，准备输出目录，启动多线程处理。
    """
    for clr in ['g', 'r', 'i']:
        check_dir(IMG_DIR_2+f"{clr}_band/", clean=True)

    parms = {
        'image_dir':        IMG_DIR,
        'image_out_dir':    IMG_DIR_2,
        'isotab_dir':       LIMIT_DEPTH_ISOTAB      
    }

    run_multi(run_single, TABLE_PATH, parms, ncpu=120)

if __name__ == "__main__":

    run_all()