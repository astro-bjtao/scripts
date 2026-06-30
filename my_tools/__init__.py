#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
my_tools — StellarHalo Process2 通用工具集

子模块
------
* ``_io``       — 目录操作、多进程调度
* ``_image``    — 图像几何变换、形态学、显示拉伸
* ``_plot``     — 比例尺、目视检查图、饱和检查三联图
* ``_cosmo``    — 宇宙学距离换算
* ``_bkg``      — 背景测量
* ``_fitting``  — Ellipse 拟合、矩估计、等强度线平滑

所有公开函数均从顶层 ``my_tools`` 直接导入，保持向后兼容。
"""

# --- io ---
from my_tools._io import check_dir, run_multi

# --- image ---
from my_tools._image import (
    get_central, cut_image, bin_image_2x2, cut_bin_image,
    smooth_img, get_circle, top_hat_kernel, extend_mask,
    get_ellipse, reg_to_mask, Lognorm2, Band2RGB, make_masked_image,
)

# --- cosmo ---
from my_tools._cosmo import (
    COSMO, kpc_per_arcsec, kpc_per_pixels, pixels_per_kpc,
    convert_arcsec2kpc, Pixs_1kpc,
)

# --- bkg ---
from my_tools._bkg import Intens2SB, MeasureBkg

# --- fitting ---
from my_tools._fitting import (
    convert_wcs, Ellipse_free, In_Ellipse,
    iso_to_table, reshape_isotable,
    get_initsma, try_fit,
    moments_estimate, smooth_isotable,
)

# --- plot ---
from my_tools._plot import plot_scalebar, plot_eyeball, plot_eyeball_saturate
