#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
根据圆形面亮度轮廓，找出背景区域.
"""

from scipy.optimize import curve_fit
from astropy.table import Table, Column
from astropy.io import fits
import numpy as np
import shutil, os

# 路径配置
from config import *


def FitLinear(x, y):
    # y = a*x + b
    def func(x, a, b):
        return a*x + b

    # initial guesses for a,b,c:
    p0 = 0., 0.
    a, b = curve_fit(func, x, y, p0)[0]
    return a, b


def FindBkg(isotab):
    # available
    max_intens = 0.01123 # 28 mag arcsec^-2
    avai   = (~np.isnan(isotab['intens'].data)) & (isotab['intens'].data < max_intens)
    isotab = isotab[avai]
    ### find bkg
    sma    = isotab['sma'].data
    intens = isotab['intens'].data
    length = len(isotab)
    index  = 0
    list_a     = []
    list_index = []
    # index bkg
    while True:
        # edge
        if (index + 6) > length:
            index -= 1 # !
            break
        # fitlinear from 6 annuli
        xdata = sma[index:(index+6)]
        ydata = intens[index:(index+6)]
        a, b = FitLinear(xdata, ydata)
        list_a.append(a)
        list_index.append(index)
        index += 1
    # esitmate bkg, 6 annuli, 24 regions.
    for i in range(len(list_a)-2):
        index_bkg = i
        if list_a[i] > 0:
            index_bkg = list_index[i]
            break
        if (list_a[i] > list_a[i+1]) & (list_a[i] > list_a[i+2]):
            index_bkg = list_index[i]
            break
    sma_min = sma[index_bkg]

    return sma_min

def GetLocalBkg():

    ## paths
    dir_iso    = LIMIT_DEPTH_ISOTAB

    # table
    table = Table.read(TABLE_PATH)
    # indexs
    for i in range(len(table)):
        # information
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]
        ### get bkgrmin
        # path
        suffix_clr = f"a_band/{label}_{index}.fits"
        path_iso = os.path.join(dir_iso, suffix_clr)
        isotab = Table.read(path_iso)
        sma_min = FindBkg(isotab)
        # save bkgr
        table['bkg_rmin'][i] = sma_min
        table['bkg_rmax'][i] = sma_min * 1.1**6
    # save table
    table.write(TABLE_PATH, overwrite=True)

if __name__ == "__main__":

    GetLocalBkg()