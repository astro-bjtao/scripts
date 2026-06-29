#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
画出拟合的等轮廓线。
"""

from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy.table import Table
from astropy.io import fits

import numpy as np
import shutil, os

from astropy.visualization import ImageNormalize
from astropy.visualization.stretch import LogStretch, LinearStretch
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
#plt.style.use('classic')

### set rcParams
plt.rcParams['font.sans-serif'] = "Times New Roman"
plt.rcParams['font.serif'] = "Times New Roman"
config = {
    "font.family":'serif',
    'font.size': 25,
    'mathtext.fontset':'stix',
}
plt.rcParams.update(config)
plt.rcParams['axes.labelsize'] = 25
plt.rcParams['figure.figsize'] = (8, 8)

# 路径配置
from config import *
# 通用工具
from my_tools import *


def Show_deepimg_1re(table, ind, dir_img, dir_mask, isotab, ax, mag_norm=20):
    # information
    label = str(table['survey'].data[ind], encoding='utf-8')
    index = table['index'].data[ind]
    redshift  = table['redshift'].data[ind]
    # path
    path_img_g = os.path.join(dir_img, f"g_band/{label}_{index}.fits")
    path_img_r = os.path.join(dir_img, f"r_band/{label}_{index}.fits")
    path_img_i = os.path.join(dir_img, f"i_band/{label}_{index}.fits")
    path_img_a = os.path.join(dir_img, f"a_band/{label}_{index}.fits")
    path_mask  = os.path.join(dir_mask,  f"{label}_{index}.fits")
    # read
    image_g = fits.getdata(path_img_g)
    image_r = fits.getdata(path_img_r)
    image_i = fits.getdata(path_img_i)
    image_a = fits.getdata(path_img_a)
    mask    = fits.getdata(path_mask)
    # cut and bin image
    cenx = int(isotab['x0'].data[0])
    ceny = int(isotab['y0'].data[0])
    sw = int(np.min([cenx, ceny])/2)
    image_g = cut_bin_image(image_g, cenx, ceny, sw)
    image_r = cut_bin_image(image_r, cenx, ceny, sw)
    image_i = cut_bin_image(image_i, cenx, ceny, sw)
    image_a = cut_bin_image(image_a, cenx, ceny, sw)
    mask    = cut_bin_image(mask,    cenx, ceny, sw)
    # smooth image wide r
    #image_a = smooth_img(image_a, std=1)
    # LogStretch for rgb
    norm_rgb = 10**((mag_norm - 27 - 5*np.log10(0.168))/(-2.5))
    image_g = Lognorm2(image_g, norm_rgb, a=1000)
    image_r = Lognorm2(image_r, norm_rgb, a=1000)
    image_i = Lognorm2(image_i, norm_rgb, a=1000)
    # LogStretch for wide r
    norm_wider = 10**((25 - 27 - 5*np.log10(0.168))/(-2.5))
    image_a = Lognorm2(image_a, norm_wider, a=1000)
    image_a_reverse = 1 - image_a
    # combine rgb and all
    flag_lowsnr = image_a < 1 # mu > 25 magqasec
    image_g[flag_lowsnr] = image_a_reverse[flag_lowsnr]
    image_r[flag_lowsnr] = image_a_reverse[flag_lowsnr]
    image_i[flag_lowsnr] = image_a_reverse[flag_lowsnr]
    # mask image
    image_g[mask>0] = 1
    image_r[mask>0] = 1
    image_i[mask>0] = 1
    # rgb image
    rgb_image = Band2RGB(image_i, image_r, image_g)
    
    # show image
    ax.imshow(rgb_image, origin='lower', interpolation='none')
    ax.axis('off')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    # scale bar
    tot_pix = 2*sw
    plot_scalebar(ax, tot_pix, redshift)

def Show_deepimg_2re(table, ind, dir_img, dir_mask, isotab, ax, mag_norm=20):
    # information
    label = str(table['survey'].data[ind], encoding='utf-8')
    index = table['index'].data[ind]
    redshift  = table['redshift'].data[ind]
    # path
    path_img_g = os.path.join(dir_img, f"g_band/{label}_{index}.fits")
    path_img_r = os.path.join(dir_img, f"r_band/{label}_{index}.fits")
    path_img_i = os.path.join(dir_img, f"i_band/{label}_{index}.fits")
    path_img_a = os.path.join(dir_img, f"a_band/{label}_{index}.fits")
    path_mask  = os.path.join(dir_mask,  f"{label}_{index}.fits")
    # read
    image_g = fits.getdata(path_img_g)
    image_r = fits.getdata(path_img_r)
    image_i = fits.getdata(path_img_i)
    image_a = fits.getdata(path_img_a)
    mask    = fits.getdata(path_mask)
    # cut and bin image
    cenx = int(isotab['x0'].data[0])
    ceny = int(isotab['y0'].data[0])
    sw = int(np.min([cenx, ceny])/4)
    image_g = cut_image(image_g, cenx, ceny, sw)
    image_r = cut_image(image_r, cenx, ceny, sw)
    image_i = cut_image(image_i, cenx, ceny, sw)
    image_a = cut_image(image_a, cenx, ceny, sw)
    mask    = cut_image(mask,    cenx, ceny, sw)
    # smooth image wide r
    #image_a = smooth_img(image_a, std=1)
    # LogStretch for rgb
    norm_rgb = 10**((mag_norm - 27 - 5*np.log10(0.168))/(-2.5))
    image_g = Lognorm2(image_g, norm_rgb, a=1000)
    image_r = Lognorm2(image_r, norm_rgb, a=1000)
    image_i = Lognorm2(image_i, norm_rgb, a=1000)
    # LogStretch for wide r
    norm_wider = 10**((25 - 27 - 5*np.log10(0.168))/(-2.5))
    image_a = Lognorm2(image_a, norm_wider, a=1000)
    image_a_reverse = 1 - image_a
    # combine rgb and all
    flag_lowsnr = image_a < 1 # mu > 25 magqasec
    image_g[flag_lowsnr] = image_a_reverse[flag_lowsnr]
    image_r[flag_lowsnr] = image_a_reverse[flag_lowsnr]
    image_i[flag_lowsnr] = image_a_reverse[flag_lowsnr]
    # mask image
    image_g[mask>0] = 1
    image_r[mask>0] = 1
    image_i[mask>0] = 1
    # rgb image
    rgb_image = Band2RGB(image_i, image_r, image_g)
    
    # show image
    ax.imshow(rgb_image, origin='lower', interpolation='none')
    ax.axis('off')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    # scale bar
    tot_pix = 2*sw
    plot_scalebar(ax, tot_pix, redshift)

def plot_patch(ax, x0, y0, sma, eps, pa_deg, color='r', lw=1, fill=False):
    # 计算短轴和角度转换
    b = sma * (1 - eps)
    width = sma # binning, normal: 2 * sma
    height = b  # normal: 2 * b
    
    # 创建椭圆对象
    ell = patches.Ellipse((x0, y0), width, height, angle=pa_deg, 
                          edgecolor=color, facecolor='none' if not fill else color,
                          lw=lw)
    ax.add_patch(ell)

def Show_smb05_1re(table, ind, isotab, ax):
    # information
    r27  = table['proc1_r27'].data[ind]
    # parms
    cenx = int(isotab['x0'].data[0])
    ceny = int(isotab['y0'].data[0])
    sw = int(np.min([cenx, ceny])/2)
    # plot isotab
    #minsma = 1.0 / pix2rekpc(redshift, re_kpc) # 1 mean re in pixel
    minsma = r27 / 2
    isotab = isotab[isotab['sma'] > minsma]
    # stop code != 4
    flag_stop = (isotab['stop_code'] != 4)
    isotab = isotab[flag_stop]
    # isotab, ind_iso
    for ind_iso in range(len(isotab)):
        sma = isotab['sma'].data[ind_iso]
        eps = isotab['ellipticity'].data[ind_iso]
        pa  = isotab['pa'].data[ind_iso]
        plot_patch(ax, sw/2, sw/2, sma, eps, pa, color='r', lw=1, fill=False)

def Show_smb05_2re(table, ind, isotab, ax):
    # information
    r27  = table['proc1_r27'].data[ind]
    # parms
    cenx = int(isotab['x0'].data[0])
    ceny = int(isotab['y0'].data[0])
    sw = int(np.min([cenx, ceny])/4)
    # plot isotab
    #minsma = 0.1 / pix2rekpc(redshift, re_kpc) # 1 mean re in pixel
    minsma = r27 / 20
    isotab = isotab[isotab['sma'] > minsma]
    maxsma = r27
    isotab = isotab[isotab['sma'] < maxsma]
    # stop code != 4
    #flag_stop = (isotab['stop_code'] != 4)
    #isotab = isotab[flag_stop]
    # isotab, ind_iso
    for ind_iso in range(len(isotab)):
        sma = isotab['sma'].data[ind_iso]
        eps = isotab['ellipticity'].data[ind_iso]
        pa  = isotab['pa'].data[ind_iso]
        plot_patch(ax, sw, sw, sma*2, eps, pa, color='r', lw=1, fill=False)

def Show_stop_4(table, ind, isotab, ax):
    # parms
    cenx = int(isotab['x0'].data[0])
    ceny = int(isotab['y0'].data[0])
    sw = int(np.min([cenx, ceny])/2)
    # stop code != 4
    flag_stop = (isotab['stop_code'] == 4) & (isotab['sma'] < sw/1.1)
    isotab = isotab[flag_stop]
    if len(isotab) < 1:
        return
    # isotab, ind_iso
    for ind_iso in range(len(isotab)):
        sma = isotab['sma'].data[ind_iso]
        eps = isotab['ellipticity'].data[ind_iso]
        pa  = isotab['pa'].data[ind_iso]
        plot_patch(ax, sw/2, sw/2, sma, eps, pa, color='gold', lw=1, fill=False)

def Show_info(table, ind, ax):
    # information
    label = str(table['survey'].data[ind], encoding='utf-8')
    index = table['index'].data[ind]
    totmass   = table['logmass_survey'].data[ind]

    # information text
    bbox = {"facecolor": "white", "alpha": 1.0}
    styles = {"size": 25, "color": "black", "bbox": bbox}
    t = f"{label.upper()}: {index}\nmass: {totmass:.2f}"
    ax.text(0.03,0.97, t, transform=ax.transAxes, 
            ha='left', va="top", linespacing=1.1, **styles)

### geometry

def show_isophotes(table, ind, dir_img, dir_mask, isotab, path_show):
    # plot img
    fig = plt.figure(figsize=(16, 8), dpi=300)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1], figure=fig)
    
    plt.subplots_adjust(top=1, bottom=0,right=1,left=0, hspace=0, wspace=0)
    plt.margins(0, 0)
    plt.xticks([])
    plt.yticks([])

    # Create main image area
    ax0 = fig.add_subplot(gs[0]) 
    ax1 = fig.add_subplot(gs[1]) 

    Show_deepimg_1re(table, ind, dir_img, dir_mask, isotab, ax0, mag_norm=20)
    Show_deepimg_2re(table, ind, dir_img, dir_mask, isotab, ax1, mag_norm=20)

    Show_smb05_1re(table, ind, isotab, ax0)
    Show_smb05_2re(table, ind, isotab, ax1)

    #Show_stop_4(table, ind, isotab, ax0)

    Show_info(table, ind, ax0)
    
    fig.savefig(path_show, bbox_inches="tight", pad_inches=0)
    plt.close()


### fitting and show result

def run_single(table, parms):

    # parms
    dir_img     = parms['img_dir']
    dir_mask    = parms['mask_dir']
    dir_iso     = parms['iso_dir']
    dir_eyeball = parms['eyeball_dir']

    for ind in range(len(table)):
        # information
        label = str(table['survey'].data[ind], encoding='utf-8')
        index = table['index'].data[ind]
        # paths
        path_iso     = os.path.join(dir_iso,     f"{label}_{index}.fits")
        path_eyeball = os.path.join(dir_eyeball, f"{label}_{index}.png")
        try:
            iso_table = Table.read(path_iso)
            show_isophotes(table, ind, dir_img, dir_mask, iso_table, path_eyeball)
        except:
            print(f"{label}_{index} failed to fit isophotes")


def run_all():
    """
    主入口：设置路径，准备输出目录，启动多线程处理。
    """
    check_dir(PROPS_ORIGINAL_EYEBALL_FITTING, clean=True)

    parms = {
        'img_dir':      IMG_DIR_2,
        'mask_dir':     TOTAL_MASK_TOTAL,
        'iso_dir':      PROPS_ORIGINAL_ISOTAB,
        'eyeball_dir':  PROPS_ORIGINAL_EYEBALL_FITTING      
    }

    run_multi(run_single, TABLE_PATH, parms, ncpu=120)

if __name__ == "__main__":

    run_all()