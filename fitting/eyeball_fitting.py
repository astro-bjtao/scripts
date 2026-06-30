#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
画出拟合的等轮廓线。
"""

from astropy.table import Table
from astropy.io import fits

import numpy as np
import os

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches

### set rcParams
plt.rcParams['font.sans-serif'] = "Times New Roman"
plt.rcParams['font.serif'] = "Times New Roman"
config = {
    "font.family": 'serif',
    'font.size': 25,
    'mathtext.fontset': 'stix',
}
plt.rcParams.update(config)
plt.rcParams['axes.labelsize'] = 25
plt.rcParams['figure.figsize'] = (8, 8)

# 路径配置
from config import *
# 通用工具
from my_tools import *


# ============================================================
#  画图组件
# ============================================================

def load_galaxy_images(table, ind, dir_img):
    """读取 g/r/i/a 波段图像和 mask，返回 dict"""
    label = str(table['survey'].data[ind], encoding='utf-8')
    index = table['index'].data[ind]
    images = {}
    for band in ['g', 'r', 'i', 'a']:
        path = os.path.join(dir_img, f"{band}_band/{label}_{index}.fits")
        images[band] = fits.getdata(path)
    return label, index, images


def show_deepimg(images, mask, cenx, ceny, redshift, ax, zoom=1, mag_norm=20):
    """
    绘制深场 RGB 图 + 比例尺。

    Parameters
    ----------
    zoom : int
        1 = 1Re 视场（bin 2×2，宽视场低表面亮度）
        2 = 2Re 视场（原始分辨率，中心区域）
    """
    if zoom == 1:
        sw = int(np.min([cenx, ceny]) / 2)
        g = cut_bin_image(images['g'], cenx, ceny, sw)
        r = cut_bin_image(images['r'], cenx, ceny, sw)
        i = cut_bin_image(images['i'], cenx, ceny, sw)
        a = cut_bin_image(images['a'], cenx, ceny, sw)
        m = cut_bin_image(mask,         cenx, ceny, sw)
    else:  # zoom == 2
        sw = int(np.min([cenx, ceny]) / 4)
        g = cut_image(images['g'], cenx, ceny, sw)
        r = cut_image(images['r'], cenx, ceny, sw)
        i = cut_image(images['i'], cenx, ceny, sw)
        a = cut_image(images['a'], cenx, ceny, sw)
        m = cut_image(mask,         cenx, ceny, sw)

    # 对数拉伸
    norm_rgb = 10 ** ((mag_norm - 27 - 5 * np.log10(0.168)) / (-2.5))
    norm_a   = 10 ** ((25 - 27 - 5 * np.log10(0.168)) / (-2.5))

    g = Lognorm2(g, norm_rgb, a=1000)
    r = Lognorm2(r, norm_rgb, a=1000)
    i = Lognorm2(i, norm_rgb, a=1000)
    a = Lognorm2(a, norm_a,   a=1000)
    a_rev = 1 - a

    # 低表面亮度区域用 a 波段灰度替代
    flag_lowsnr = a < 1  # mu > 25 mag/arcsec²
    g[flag_lowsnr] = a_rev[flag_lowsnr]
    r[flag_lowsnr] = a_rev[flag_lowsnr]
    i[flag_lowsnr] = a_rev[flag_lowsnr]

    # 叠加掩模（白色）
    g[m > 0] = 1
    r[m > 0] = 1
    i[m > 0] = 1

    rgb = Band2RGB(i, r, g)

    ax.imshow(rgb, origin='lower', interpolation='none')
    ax.axis('off')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    # 比例尺
    plot_scalebar(ax, 2 * sw, redshift)


def show_isophote_patches(table, ind, isotab, ax, zoom=1):
    """
    在图上叠加等强度椭圆。

    Parameters
    ----------
    zoom : int
        1 = 1Re（bin 后坐标系，sma 即像素坐标）
        2 = 2Re（原始坐标系，sma×2 以匹配 bin 后的视觉效果）
    """
    r27 = table['proc1_r27'].data[ind]

    if zoom == 1:
        sw = int(np.min([isotab['x0'].data[0], isotab['y0'].data[0]]) / 2)
        minsma = r27 / 2
        maxsma = None
        sma_scale = 1.0
        patch_cx, patch_cy = sw / 2, sw / 2
        # 只保留 stop_code != 4 的椭圆
        isotab = isotab[isotab['stop_code'] != 4]
    else:  # zoom == 2
        sw = int(np.min([isotab['x0'].data[0], isotab['y0'].data[0]]) / 4)
        minsma = r27 / 20
        maxsma = r27
        sma_scale = 2.0
        patch_cx, patch_cy = sw, sw

    isotab = isotab[isotab['sma'] > minsma]
    if maxsma is not None:
        isotab = isotab[isotab['sma'] < maxsma]

    for row in isotab:
        plot_patch(ax, patch_cx, patch_cy,
                   row['sma'] * sma_scale, row['ellipticity'], row['pa'],
                   color='r', lw=1, fill=False)


def show_stop_4(table, ind, isotab, ax):
    """叠加 stop_code==4 的椭圆（金色），用于诊断拟合异常终止的位置"""
    cenx = int(isotab['x0'].data[0])
    ceny = int(isotab['y0'].data[0])
    sw = int(np.min([cenx, ceny]) / 2)

    flag = (isotab['stop_code'] == 4) & (isotab['sma'] < sw / 1.1)
    isotab = isotab[flag]
    if len(isotab) < 1:
        return

    for row in isotab:
        plot_patch(ax, sw / 2, sw / 2,
                   row['sma'], row['ellipticity'], row['pa'],
                   color='gold', lw=1, fill=False)


def show_info(table, ind, ax):
    """左上角标注星系 ID 和质量"""
    label = str(table['survey'].data[ind], encoding='utf-8')
    index = table['index'].data[ind]
    totmass = table['logmass_survey'].data[ind]

    bbox = {"facecolor": "white", "alpha": 1.0}
    styles = {"size": 25, "color": "black", "bbox": bbox}
    t = f"{label.upper()}: {index}\nmass: {totmass:.2f}"
    ax.text(0.03, 0.97, t, transform=ax.transAxes,
            ha='left', va="top", linespacing=1.1, **styles)


def plot_patch(ax, x0, y0, sma, eps, pa_deg, color='r', lw=1, fill=False):
    """在 ax 上绘制一个椭圆 patch"""
    b = sma * (1 - eps)
    ell = patches.Ellipse((x0, y0), sma, b, angle=pa_deg,
                          edgecolor=color,
                          facecolor='none' if not fill else color,
                          lw=lw)
    ax.add_patch(ell)


# ============================================================
#  主图
# ============================================================

def show_isophotes(table, ind, dir_img, dir_mask, isotab, path_show):
    """绘制 1Re + 2Re 对比图，叠加等强度线"""
    # 提前加载所有图像数据（避免 1re/2re 重复读取）
    label, index, images = load_galaxy_images(table, ind, dir_img)
    mask = fits.getdata(os.path.join(dir_mask, f"{label}_{index}.fits"))
    redshift = table['redshift'].data[ind]
    cenx = int(isotab['x0'].data[0])
    ceny = int(isotab['y0'].data[0])

    fig = plt.figure(figsize=(16, 8), dpi=300)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1], figure=fig)

    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    plt.margins(0, 0)
    plt.xticks([])
    plt.yticks([])

    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])

    # 1Re 视场（bin 2×2，宽视场）
    show_deepimg(images, mask, cenx, ceny, redshift, ax0, zoom=1, mag_norm=20)
    show_isophote_patches(table, ind, isotab, ax0, zoom=1)

    # 2Re 视场（原始分辨率，中心区域）
    show_deepimg(images, mask, cenx, ceny, redshift, ax1, zoom=2, mag_norm=20)
    show_isophote_patches(table, ind, isotab, ax1, zoom=2)

    show_info(table, ind, ax0)

    fig.savefig(path_show, bbox_inches="tight", pad_inches=0)
    plt.close()


# ============================================================
#  多进程调度
# ============================================================

def run_single(table, parms):
    dir_img     = parms['img_dir']
    dir_mask    = parms['mask_dir']
    dir_iso     = parms['iso_dir']
    dir_eyeball = parms['eyeball_dir']

    for ind in range(len(table)):
        label = str(table['survey'].data[ind], encoding='utf-8')
        index = table['index'].data[ind]
        path_iso     = os.path.join(dir_iso,     f"{label}_{index}.fits")
        path_eyeball = os.path.join(dir_eyeball, f"{label}_{index}.png")

        try:
            iso_table = Table.read(path_iso)
            show_isophotes(table, ind, dir_img, dir_mask, iso_table, path_eyeball)
        except Exception:
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
        'eyeball_dir':  PROPS_ORIGINAL_EYEBALL_FITTING,
    }

    run_multi(run_single, TABLE_PATH, parms, ncpu=120)


if __name__ == "__main__":
    run_all()
