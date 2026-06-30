#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
背景测量工具。
"""

import numpy as np


def Intens2SB(intens, A=27, pixscale=0.168):
    """
    流量强度 → 面亮度（mag/arcsec²）转换。

    SB = -2.5·log10(intens) + A + 5·log10(pixscale)
    """
    return -2.5 * np.log10(intens) + A + 5 * np.log10(pixscale)


def MeasureBkg(isotab, sma_min):
    """
    在背景环形区域中测量 rms、高低频噪声和局部背景值。

    使用 sma_min ~ sma_max（1.1^6 倍范围）内的 6 个椭圆环，
    每环 4 个象限共 24 个数据点估计背景统计量。

    Parameters
    ----------
    isotab : Table
        四象限等强度线表（含 intens_{qx}, intens_err_{qx}, rms_{qx} 列）。
    sma_min : float
        背景区域起始半长轴。

    Returns
    -------
    rms : float
        中位数 RMS（像素值单位）。
    err_high : float
        高频噪声（quadrant 间中位数误差）。
    err_low : float
        低频噪声（quadrant 间通量 scatter 扣除高频分量后）。
    err_tot : float
        总误差 = sqrt(high² + low²)。
    local_bkg : float
        局部背景值（quadrant 中位数通量）。
    """
    sma_max = sma_min * 1.1 ** 6
    semi_step = 1.1 ** 0.5
    flag = (sma_min / semi_step < isotab['sma'].data) & \
           (isotab['sma'].data < sma_max * semi_step)
    iso = isotab[flag]

    # 收集 24 个象限的数据
    intens_list, interr_list, rms_list = [], [], []
    for qx in range(4):
        for i in range(6):
            intens_list.append(iso[f'intens_{qx}'].data[i])
            interr_list.append(iso[f'intens_err_{qx}'].data[i])
            rms_list.append(iso[f'rms_{qx}'].data[i])

    rms = np.nanmedian(rms_list)
    err_high = np.nanmedian(interr_list)
    var_high = err_high ** 2

    std_sky = np.nanstd(intens_list)
    var_low = max(0.0, std_sky ** 2 - var_high)
    err_low = np.sqrt(var_low)
    err_tot = np.sqrt(var_high + var_low)

    local_bkg = np.nanmedian(intens_list)

    return rms, err_high, err_low, err_tot, local_bkg
