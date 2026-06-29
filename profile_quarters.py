#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
根据减过亮星的图像, mask_total, 来获取背景区域面亮度分布.
"""

from astropy.table import Table, Column
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import shutil, os
import math

# 路径配置
from config import *
# 通用工具
from my_tools import check_dir, run_multi

### Profile Manual
def DistCircle(shape, cenx, ceny):

    temp = np.indices(shape, dtype=np.float32)
    xtp = temp[1] - cenx
    ytp = temp[0] - ceny
    im = np.hypot(xtp, ytp)
    im_pa = np.arctan2(xtp,ytp)*180/math.pi+180
    return im, im_pa

def RegionCircle(shape, cenx, ceny, minr, maxr):

    im, im_pa = DistCircle(shape, cenx, ceny)
    return  (minr <= im) & (im < maxr)

def RegionQuarter(shape, cenx, ceny, minr, maxr, qx):
    # qx: 0, 1, 2, 3
    im, im_pa = DistCircle(shape, cenx, ceny)
    return  (minr <= im) & (im < maxr) & (qx*90 <= im_pa) & (im_pa < (qx+1)*90)

def intens_region(data, region, sclip=5.0, nclip=5):
    
    ndata_a = np.count_nonzero((~np.isnan(data)) & region) # avaiable ndata
    ndata_t = region.sum() # total ndata
    data = data[region] # cut data
    # return if ndata_a is less than 5
    if ndata_a < 5:
        intens=np.nan
        stddev=np.nan
        int_err=np.nan
        return intens, stddev, int_err, ndata_a, ndata_t
    # return if ndata_a is less than 0.1
    if (ndata_a/ndata_t < 0.1):
        intens=np.nan
        stddev=np.nan
        int_err=np.nan
        return intens, stddev, int_err, ndata_a, ndata_t
    # mean or sigma clipping mean
    mean, median, stddev = sigma_clipped_stats(data, mask=None, mask_value=None, 
                                               sigma=sclip, sigma_lower=None, sigma_upper=None, maxiters=nclip, 
                                               cenfunc='median', stdfunc='std', std_ddof=0, axis=None, grow=False)
    intens = median
    int_err = stddev / np.sqrt(ndata_a) # error of intensity

    return intens, stddev, int_err, ndata_a, ndata_t

def AppendProperties(data, region, list_prop):

    list_intens, list_int_err, list_rms, list_ndata_a, list_ndata_t = list_prop
    intens, rms, int_err, ndata_a, ndata_t = intens_region(data, region)
    list_intens.append(intens)
    list_int_err.append(int_err)
    list_rms.append(rms)
    list_ndata_a.append(ndata_a)
    list_ndata_t.append(ndata_t)

def List2Column(tab, list_prop):

    list_intens, list_int_err, list_rms, list_ndata_a, list_ndata_t = list_prop
    # array
    intens = np.array(list_intens)
    int_err = np.array(list_int_err)
    rms = np.array(list_rms)
    ndata_a = np.array(list_ndata_a)
    ndata_t = np.array(list_ndata_t)
    # column
    intens = Column(data=intens, name='intens', dtype=np.float32)
    int_err = Column(data=int_err, name='intens_err', dtype=np.float32)
    rms = Column(data=rms, name='rms', dtype=np.float32)
    ndata_a = Column(data=ndata_a, name='ndata_a', dtype=np.float32)
    ndata_t = Column(data=ndata_t, name='ndata_t', dtype=np.float32)
    # add columns
    tab.add_columns([intens, int_err, rms, ndata_a, ndata_t])

    return tab

def List2ColumnQX(tab, list_prop, qx):

    list_intens, list_int_err, list_rms, list_ndata_a, list_ndata_t = list_prop
    # array
    intens = np.array(list_intens)
    int_err = np.array(list_int_err)
    rms = np.array(list_rms)
    ndata_a = np.array(list_ndata_a)
    ndata_t = np.array(list_ndata_t)
    # column
    intens = Column(data=intens, name=f'intens_{qx}', dtype=np.float32)
    int_err = Column(data=int_err, name=f'intens_err_{qx}', dtype=np.float32)
    rms = Column(data=rms, name=f'rms_{qx}', dtype=np.float32)
    ndata_a = Column(data=ndata_a, name=f'ndata_a_{qx}', dtype=np.float32)
    ndata_t = Column(data=ndata_t, name=f'ndata_t_{qx}', dtype=np.float32)
    # add columns
    tab.add_columns([intens, int_err, rms, ndata_a, ndata_t])

    return tab


def ProfileQuarters(data, xc, yc, max_sma):
    """
    data: masked image, in which the mask is nan
    xc, yc: center
    max_sma: maximum semi-major axis
    """
    ### initial
    list_sma = []
    # total
    list_intens = []
    list_int_err = []
    list_rms = []
    list_ndata_a = []
    list_ndata_t = []
    list_total = [list_intens, list_int_err, list_rms, list_ndata_a, list_ndata_t]
    # q0
    list0_intens = []
    list0_int_err = []
    list0_rms = []
    list0_ndata_a = []
    list0_ndata_t = []
    list_q0 = [list0_intens, list0_int_err, list0_rms, list0_ndata_a, list0_ndata_t]
    # q1
    list1_intens = []
    list1_int_err = []
    list1_rms = []
    list1_ndata_a = []
    list1_ndata_t = []
    list_q1 = [list1_intens, list1_int_err, list1_rms, list1_ndata_a, list1_ndata_t]
    # q2
    list2_intens = []
    list2_int_err = []
    list2_rms = []
    list2_ndata_a = []
    list2_ndata_t = []
    list_q2 = [list2_intens, list2_int_err, list2_rms, list2_ndata_a, list2_ndata_t]
    # q3
    list3_intens = []
    list3_int_err = []
    list3_rms = []
    list3_ndata_a = []
    list3_ndata_t = []
    list_q3 = [list3_intens, list3_int_err, list3_rms, list3_ndata_a, list3_ndata_t]
    # list tuple qx
    list_qx = [list_q0, list_q1, list_q2, list_q3]
    ### linear
    sma = 0.5
    semi_step = 0.5
    while True:
        min_r = sma - semi_step
        max_r = sma + semi_step
        # sma
        list_sma.append(sma)
        # total
        region = RegionCircle(data.shape, xc, yc, min_r, max_r)
        AppendProperties(data, region, list_total)
        del region
        # quarters
        for qx in range(4):
            region = RegionQuarter(data.shape, xc, yc, min_r, max_r, qx)
            AppendProperties(data, region, list_qx[qx])
            del region
        # next
        sma = sma + 1
        if sma >= 10:
            break
    ### no-linear
    sma = 10.55688837
    semi_step = (1.1)**0.5
    while True:
        min_r = sma / semi_step
        max_r = sma * semi_step
        # sma
        list_sma.append(sma)
        # total
        region = RegionCircle(data.shape, xc, yc, min_r, max_r)
        AppendProperties(data, region, list_total)
        del region
        # quarters
        for qx in range(4):
            region = RegionQuarter(data.shape, xc, yc, min_r, max_r, qx)
            AppendProperties(data, region, list_qx[qx])
            del region
        # next
        sma = sma * 1.1
        if sma >= max_sma:
            break
    ### Table
    tab = Table()
    # sma
    sma = np.array(list_sma)
    sma = Column(data=sma, name='sma', dtype=np.float32)
    tab.add_columns([sma,])
    # total
    tab = List2Column(tab, list_total)
    # quarters
    for qx in range(4):
        tab = List2ColumnQX(tab, list_qx[qx], qx)
    
    return tab


### geometry

def convert_wcs(ra, dec, hdr):

    sky = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    w = WCS(header=hdr)
    xc, yc = w.world_to_pixel(sky)
    xc = float(xc)
    yc = float(yc)

    return xc, yc

# mask thresh multiprocess

def run_single(table, parms):

    dir_img  = parms['image_dir']
    dir_mask = parms['mask_dir']
    dir_iso  = parms['isotab_dir']

    for i in range(len(table)):
        # information
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]
        ra    = table['ra'].data[i]
        dec   = table['dec'].data[i]
        # path
        for clr in ['g', 'r', 'i', 'a']:
            # paths
            suffix     = f"{label}_{index}.fits"
            suffix_clr = f"{clr}_band/{label}_{index}.fits"
            path_img   = os.path.join(dir_img,  suffix_clr)
            path_mask  = os.path.join(dir_mask, suffix)
            path_iso   = os.path.join(dir_iso,  suffix_clr)
            # image
            img, hdr = fits.getdata(path_img, header=True)
            mask = fits.getdata(path_mask)
            img[mask>0] = np.nan
            # geometry
            xc, yc = convert_wcs(ra, dec, hdr)
            sw = np.min(img.shape)/2
            isotab = ProfileQuarters(img, xc, yc, sw)
            isotab.write(path_iso)

def run_all():
    """
    主入口：设置路径，准备输出目录，启动多线程处理。
    """
    for clr in ['g', 'r', 'i', 'a']:
        check_dir(LIMIT_DEPTH_ISOTAB+f"{clr}_band/", clean=True)

    parms = {
        'image_dir':   IMG_DIR,
        'mask_dir':    TOTAL_MASK_TOTAL,
        'isotab_dir':  LIMIT_DEPTH_ISOTAB      
    }

    run_multi(run_single, TABLE_PATH, parms, ncpu=120)

if __name__ == "__main__":

    run_all()
