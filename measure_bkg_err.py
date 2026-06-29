#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
根据背景区域的1/4圈环值，获取背景误差。
"""

from astropy.table import Table, Column
from astropy.io import fits
import numpy as np
import shutil, os

# 路径配置
from config import *

### Function
# 通用工具
from my_tools import Intens2SB, MeasureBkg


def GetLocalBkg():

    ## paths
    dir_iso    = LIMIT_DEPTH_ISOTAB
    path_table = TABLE_PATH
    # table
    table = Table.read(path_table)
    # indexs
    for i in range(len(table)):
        # information
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]
        bkg_rmin = table['bkg_rmin'].data[i]
        ### get bkgrmin
        # path
        for clr in ['g', 'r', 'i']:
            # paths
            suffix_clr = f"{clr}_band/{label}_{index}.fits"
            path_iso    = os.path.join(dir_iso, suffix_clr)
            isotab = Table.read(path_iso)
            rms, err_high, err_low, err_tot, local_bkg = MeasureBkg(isotab, bkg_rmin)
            # save bkg properties of g, r, i band
            table[f'bkg_std_{clr}'][i] = rms
            table[f'bkg_hfqerr_{clr}'][i] = err_high
            table[f'bkg_lfqerr_{clr}'][i] = err_low
            table[f'bkg_toterr_{clr}'][i] = err_tot
            table[f'limit_mu_{clr}'][i] = Intens2SB(err_tot*2)
            table[f'bkg_valerr_{clr}'][i] = err_tot / np.sqrt(24)
    # save table
    table.write(path_table, overwrite=True)

if __name__ == "__main__":

    GetLocalBkg()