#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
宇宙学距离换算。
"""

from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

# Planck 2013 宇宙学参数
COSMO = FlatLambdaCDM(H0=67.3 * u.km / (u.s * u.Mpc), Tcmb0=2.725, Om0=0.315)


def kpc_per_arcsec(redshift):
    """给定红移下每角秒对应的 kpc 数（基于角直径距离）。"""
    dist_mpc = COSMO.angular_diameter_distance(redshift).value  # Mpc
    return dist_mpc * 4.848 / 1000  # kpc / arcsec


def kpc_per_pixels(redshift, pix_scale=0.168):
    """每像素对应的 kpc 数。pix_scale 单位 arcsec/pixel。"""
    return pix_scale * kpc_per_arcsec(redshift)


def pixels_per_kpc(redshift, pix_scale=0.168):
    """每 kpc 对应的像素数。"""
    return 1.0 / kpc_per_pixels(redshift, pix_scale)


def convert_arcsec2kpc(redshift):
    """角秒 → kpc 换算系数。"""
    return kpc_per_arcsec(redshift)


def Pixs_1kpc(redshift):
    """1 kpc 对应的像素数（pix_scale=0.168）。"""
    return pixels_per_kpc(redshift, pix_scale=0.168)
