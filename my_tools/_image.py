#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
图像几何变换、形态学操作、显示拉伸。
"""

import math
import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve2d
from regions import Regions


# ============================================================
#  几何操作
# ============================================================

def get_central(shape):
    """返回图像中心像素坐标 (cenx, ceny)，0-indexed。"""
    ny, nx = shape[0], shape[1]
    return int((nx - 1) / 2), int((ny - 1) / 2)


def cut_image(data, cenx, ceny, sw):
    """以 (cenx, ceny) 为中心、半径 sw 裁切正方形区域（无 binning）。"""
    return data[ceny - sw:ceny + sw + 1, cenx - sw:cenx + sw + 1]


def bin_image_2x2(data):
    """2×2 像素合并（均值），尺寸减半，提升低表面亮度信噪比。"""
    ny, nx = data.shape
    data = data[:ny // 2 * 2, :nx // 2 * 2]
    return (data[0::2, 0::2] + data[1::2, 0::2]
            + data[0::2, 1::2] + data[1::2, 1::2]) / 4.0


def cut_bin_image(data, cenx, ceny, sw):
    """裁切后 2×2 bin。"""
    return bin_image_2x2(cut_image(data, cenx, ceny, sw))


def smooth_img(data, std=1.0):
    """高斯平滑，sigma=std 像素。"""
    kernel = Gaussian2DKernel(std).array
    return convolve2d(data, kernel, mode='same', boundary='symm')


# ============================================================
#  形态学操作
# ============================================================

def get_circle(shape, cenx, ceny):
    """生成圆形距离矩阵（到 (cenx, ceny) 的像素距离）。"""
    y, x = np.indices(shape, dtype=np.float32)
    return np.hypot(x - cenx, y - ceny)


def top_hat_kernel(radius):
    """半径为 radius 的圆形 top-hat 卷积核（float32）。"""
    size = 2 * radius + 1
    kernel = get_circle((size, size), radius, radius) < (radius + 0.1)
    return kernel.astype(np.float32)


def extend_mask(mask, ext_r=3):
    """
    二值掩模膨胀操作（卷积 + 阈值）。

    Parameters
    ----------
    mask : 2D ndarray
        输入二值掩模。
    ext_r : int
        膨胀半径（像素），默认 3。

    Returns
    -------
    extended : 2D ndarray (bool)
        膨胀后的二值掩模。
    """
    kernel = top_hat_kernel(ext_r)
    smoothed = convolve2d(mask, kernel, mode='same', boundary='symm')
    return smoothed > 0


def get_ellipse(shape, cenx, ceny, b2a, posang):
    """
    归一化椭圆距离矩阵。等于 1 时位于椭圆边界上。

    Parameters
    ----------
    shape : (ny, nx)
        图像形状。
    cenx, ceny : int
        椭圆中心。
    b2a : float
        半短轴/半长轴比 (0, 1]。
    posang : float
        位置角（度），从 X+ 逆时针。
    """
    ratio = 1.0 / b2a
    rad = math.radians(posang)
    cosang, sinang = math.cos(rad), math.sin(rad)
    y, x = np.indices(shape, dtype=np.float32)
    x_diff = x - cenx
    y_diff = y - ceny
    x_rot = x_diff * cosang + y_diff * sinang
    y_rot = -x_diff * sinang + y_diff * cosang
    return np.hypot(x_rot * ratio, y_rot)


def reg_to_mask(reg_file, reference_fits, out_mask_fits):
    """
    将 DS9 区域文件 (.reg) 转为二值掩模 FITS。

    Parameters
    ----------
    reg_file : str
        DS9 区域文件路径。
    reference_fits : str
        参考 FITS 图像（提供形状和 WCS）。
    out_mask_fits : str
        输出掩模 FITS 路径。
    """
    with fits.open(reference_fits) as hdul:
        ref_data = hdul[0].data
        header = hdul[0].header

    mask = np.zeros(ref_data.shape, dtype=np.uint8)
    regions = Regions.read(reg_file, format='ds9')

    for reg in regions:
        tmp = reg.to_mask(mode='center')
        tmp = tmp.to_image(ref_data.shape)
        mask[tmp > 0] = 1

    fits.PrimaryHDU(data=mask, header=header).writeto(out_mask_fits, overwrite=True)


# ============================================================
#  显示拉伸
# ============================================================

def Lognorm2(data, vmax, a=1000):
    """
    对数拉伸到 [0, 1]。

    Parameters
    ----------
    data : ndarray
        输入数据（会被原地修改并 clamp）。
    vmax : float
        归一化上限。
    a : float
        拉伸参数，越大越接近线性，默认 1000。
    """
    data = data / vmax
    data[data < 1e-10] = 1e-10
    data[data > 1.0] = 1.0
    y = np.log(a * data + 1) / np.log(a + 1)
    y[y > 1.0] = 1.0
    y[y < 0.0] = 0.0
    return y


def Band2RGB(img_r, img_g, img_b):
    """三通道合并为 RGB 图像。各通道应在 [0, 1] 范围。"""
    assert img_r.shape == img_g.shape == img_b.shape, "Image shapes must be identical"
    rgb = np.stack([img_r, img_g, img_b], axis=2)
    rgb[rgb < 0] = 0
    rgb[rgb > 1] = 1
    return rgb


def make_masked_image(data, mask):
    """
    用掩模遮蔽科学图像。

    Parameters
    ----------
    data : ndarray
        科学图像数据。
    mask : ndarray (bool 或 0/1)
        1/True 的像素被遮蔽。

    Returns
    -------
    np.ma.MaskedArray
    """
    return np.ma.array(data, mask=mask)
