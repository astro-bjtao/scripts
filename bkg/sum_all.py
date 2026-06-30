#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将g, r, i三个波段的图像平均成a band,或者wide r band.
"""

from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.io import fits
import numpy as np

import shutil
import os

# 路径配置
from config import *
# 通用工具
from my_tools import check_dir, run_multi

# mask thresh multiprocess

def run_single(table, parms):

    # parms
    dir_img  = parms['image_dir']

    for i in range(len(table)):
        # information
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]
        # path
        path_image_g  = os.path.join(dir_img, f"g_band/{label}_{index}.fits")
        path_image_r  = os.path.join(dir_img, f"r_band/{label}_{index}.fits")
        path_image_i  = os.path.join(dir_img, f"i_band/{label}_{index}.fits")
        path_image_a  = os.path.join(dir_img, f"a_band/{label}_{index}.fits")
        # read images
        image_g, hdr = fits.getdata(path_image_g, header=True)
        image_r = fits.getdata(path_image_r)
        image_i = fits.getdata(path_image_i)
        # sum all
        image_a = (image_g + image_r + image_i) / 3.
        # save
        fits.PrimaryHDU(data=image_a, header=hdr).writeto(path_image_a)


def run_all():
    """
    主入口：设置路径，准备输出目录，启动多线程处理。
    """
    check_dir(IMG_DIR_2+f"a_band/", clean=True)

    parms = {
        'image_dir':        IMG_DIR_2,      
    }

    run_multi(run_single, TABLE_PATH, parms, ncpu=120)

if __name__ == "__main__":

    run_all()