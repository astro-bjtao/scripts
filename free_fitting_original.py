#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
对减过背景的星系图像进行测光。
"""

from astropy.table import Table, Column
from astropy.io import fits
import numpy as np
import shutil, os

# 路径配置
from config import *
# 通用工具
from my_tools import *

### fitting and show result

def run_single(table, parms):

    # parms
    # parms
    dir_img  = parms['image_dir']
    dir_mask = parms['mask_dir']
    dir_iso  = parms['isotab_dir']

    for i in range(len(table)):
        # information
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]
        ra    = table['ra'].data[i]
        dec   = table['dec'].data[i]
        sma_27   = table['proc1_r27'].data[i]
        # paths
        suffix = f"{label}_{index}.fits"
        suffix_clr = f"a_band/{label}_{index}.fits"
        path_img     = os.path.join(dir_img,   suffix_clr)
        path_mask    = os.path.join(dir_mask,  suffix)
        path_iso     = os.path.join(dir_iso,   suffix)
        # image
        img, hdr = fits.getdata(path_img, header=True)
        mask = fits.getdata(path_mask)
        masked_image = make_masked_image(img, mask)
        # geometry
        xc, yc = convert_wcs(ra, dec, hdr)
        minsma = 1.0
        maxsma = xc/2
        list_init = (300/0.168)*1.1**(-1*np.arange(60))
        delta_init = np.abs(list_init - sma_27/2.0)
        flag_init = delta_init < (np.min(delta_init) + 0.01)
        initsma = float(list_init[flag_init][0])
        try:
            iso = Ellipse_free(masked_image, 
                                xc, yc, 0.5, 0, 
                                initsma, minsma, maxsma,
                                step=0.1, fix_center=True, fix_pa=False)
            iso_tab = iso_to_table(iso)
            iso_table = reshape_isotable(iso_tab)
            iso_table.write(path_iso)
        except:
            print(f"{label}_{index} failed to fit")

def run_all():
    """
    主入口：设置路径，准备输出目录，启动多线程处理。
    """
    check_dir(PROPS_ORIGINAL_ISOTAB, clean=True)

    parms = {
        'image_dir':   IMG_DIR_2,
        'mask_dir':    TOTAL_MASK_TOTAL,
        'isotab_dir':  PROPS_ORIGINAL_ISOTAB      
    }

    run_multi(run_single, TABLE_PATH, parms, ncpu=120)
if __name__ == "__main__":

    run_all()