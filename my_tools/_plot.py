#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
目视检查图绘制：比例尺、RGB 目视图、饱和检查三联图。
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from my_tools._image import (
    get_central, cut_image, bin_image_2x2, smooth_img, Lognorm2, Band2RGB,
)
from my_tools._cosmo import kpc_per_pixels


# ============================================================
#  比例尺
# ============================================================

def plot_scalebar(ax, tot_pix, redshift):
    """
    在 Axes 左下角绘制物理比例尺。

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    tot_pix : int
        图像总宽度（像素）。
    redshift : float
    """
    kpc_per_pix = kpc_per_pixels(redshift, pix_scale=0.168)
    totwd_kpc = tot_pix * kpc_per_pix

    # 就近选取整数值
    for val in [2, 5, 10, 20, 30, 50, 100]:
        if totwd_kpc < val * 4:
            length_bar = val
            break
    else:
        length_bar = 100

    frac = length_bar / totwd_kpc
    kw = dict(transform=ax.transAxes, color='red', lw=2)

    ax.hlines(y=0.10, xmin=0.10, xmax=0.10 + frac, **kw)
    ax.vlines(x=0.10, ymin=0.09, ymax=0.11, **kw)
    ax.vlines(x=0.10 + frac, ymin=0.09, ymax=0.11, **kw)
    ax.text(0.10 + frac / 2, 0.13, f'{length_bar} kpc',
            transform=ax.transAxes, color='red', ha='center', fontsize=12)


# ============================================================
#  RGB 目视检查图
# ============================================================

def plot_eyeball(
    img_g, img_r, img_i,
    out_path,
    mask=None,
    mode='rgb',
    norm_value=3.0,
    sw_frac=1.0,
    smooth_sigma=1.0,
    dpi=300,
    label_text=None,
):
    """
    绘制星系三色目视检查图，可选叠加掩模。

    Parameters
    ----------
    img_g, img_r, img_i : 2D ndarray
        g/r/i 波段图像。
    out_path : str
        输出 PNG 路径。
    mask : 2D ndarray or None
        二值掩模（叠加在蓝通道）。
    mode : str
        'rgb' | 'rgb_mask' | 'deep' | 'deep_mask'
    norm_value : float
        对数拉伸 vmax。
    sw_frac : float
        裁切比例（1=全图）。
    smooth_sigma : float
        deep 模式高斯平滑 σ。
    dpi : int
    label_text : str or None
        左上角标注。
    """
    cenx, ceny = get_central(img_g.shape)
    sw = int(min(cenx, ceny) * sw_frac)

    g = cut_image(img_g, cenx, ceny, sw)
    r = cut_image(img_r, cenx, ceny, sw)
    i = cut_image(img_i, cenx, ceny, sw)
    m = cut_image(mask.astype(np.float32), cenx, ceny, sw) if mask is not None else None

    if mode in ('deep', 'deep_mask'):
        g = np.abs(smooth_img(bin_image_2x2(g), std=smooth_sigma))
        r = np.abs(smooth_img(bin_image_2x2(r), std=smooth_sigma))
        i = np.abs(smooth_img(bin_image_2x2(i), std=smooth_sigma))
        if m is not None:
            m = bin_image_2x2(m) > 0.5
        gray = Lognorm2((g + r + i) / 3.0, vmax=norm_value, a=1000)

        if mode == 'deep':
            rgb = np.stack([gray] * 3, axis=2)
        else:
            bg = (1 - gray) * 0.8
            b_ch = bg + (m * 0.2 if m is not None else 0)
            rgb = Band2RGB(bg, bg, b_ch)
    else:
        g = Lognorm2(g, vmax=norm_value, a=1000)
        r = Lognorm2(r, vmax=norm_value, a=1000)
        i = Lognorm2(i, vmax=norm_value, a=1000)
        if mode == 'rgb':
            rgb = Band2RGB(i, r, g)
        else:  # rgb_mask
            g, r, i = g * 0.7, r * 0.7, i * 0.7
            if m is not None:
                i = i + m * 0.3
            rgb = Band2RGB(i, r, g)

    fig = plt.figure(figsize=(7, 7), dpi=dpi)
    fig.patch.set_visible(False)
    gs = gridspec.GridSpec(1, 1, figure=fig)
    ax = fig.add_subplot(gs[0])
    ax.imshow(rgb, origin='lower', interpolation='none')
    ax.set_xticks([]); ax.set_yticks([]); ax.axis('off')
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    plt.margins(0, 0)

    if label_text:
        ax.text(0.03, 0.93, label_text, transform=ax.transAxes,
                ha='left', va='top', fontsize=15,
                bbox={"facecolor": "white", "alpha": 0.9}, color='black')

    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)


# ============================================================
#  饱和检查三联图
# ============================================================

def plot_eyeball_saturate(
    img_g, img_r, img_i,
    out_path,
    label_text=None,
    sw_frac=1.0,
    dpi=300,
):
    """
    三联图：RGB | g-r 颜色 | g-i 颜色 + 中心流量。

    Parameters
    ----------
    img_g, img_r, img_i : 2D ndarray
    out_path : str
    label_text : str or None
    sw_frac : float
    dpi : int
    """
    cenx, ceny = get_central(img_g.shape)
    sw = int(min(cenx, ceny) * sw_frac)

    g = cut_image(img_g, cenx, ceny, sw)
    r = cut_image(img_r, cenx, ceny, sw)
    i = cut_image(img_i, cenx, ceny, sw)

    eps_ = 1e-10
    color_gr = -2.5 * np.log10(np.maximum(g, eps_) / np.maximum(r, eps_))
    color_gi = -2.5 * np.log10(np.maximum(g, eps_) / np.maximum(i, eps_))

    cval_g = g[sw, sw]; cval_r = r[sw, sw]; cval_i = i[sw, sw]
    norm = max(cval_g, cval_r, cval_i)
    rgb = Band2RGB(
        Lognorm2(i, norm, a=1000),
        Lognorm2(r, norm, a=1000),
        Lognorm2(g, norm, a=1000),
    )

    fig = plt.figure(figsize=(21, 7), dpi=dpi)
    fig.patch.set_visible(False)
    w = 1.0 / 3.0
    ax0 = fig.add_axes([0.0, 0.0, w, 1.0])
    ax1 = fig.add_axes([w, 0.0, w, 1.0])
    ax2 = fig.add_axes([2 * w, 0.0, w, 1.0])

    bbox = {"facecolor": "white", "alpha": 1.0}
    sty = {"size": 15, "color": "black", "bbox": bbox}

    ax0.imshow(rgb, origin='lower', interpolation='none', aspect='auto')
    ax0.set_xticks([]); ax0.set_yticks([]); ax0.axis('off')
    if label_text:
        ax0.text(0.03, 0.93, label_text, transform=ax0.transAxes,
                 ha='left', **sty)

    ax1.imshow(color_gr, cmap='cool', vmax=1.5, vmin=-0.5,
               origin='lower', interpolation='none', aspect='auto')
    ax1.set_xticks([]); ax1.set_yticks([]); ax1.axis('off')
    ax1.text(0.05, 0.93, "g - r color image", transform=ax1.transAxes,
             ha='left', **sty)

    ax2.imshow(color_gi, cmap='cool', vmax=2, vmin=0,
               origin='lower', interpolation='none', aspect='auto')
    ax2.set_xticks([]); ax2.set_yticks([]); ax2.axis('off')
    ax2.text(0.05, 0.93, "g - i color image", transform=ax2.transAxes,
             ha='left', **sty)
    ax2.text(0.05, 0.05,
             f"central intensity\ng: {cval_g:.1f}\nr: {cval_r:.1f}\ni: {cval_i:.1f}",
             transform=ax2.transAxes, ha='left', **sty)

    fig.savefig(out_path, dpi=dpi, pad_inches=0.0)
    plt.close(fig)
