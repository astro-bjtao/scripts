#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
计算曲线增长 (Curve of Growth) 表。

tflux_e / tflux_max → COG 值。

依赖表列: deconv_flux_{clr}, deconv_fluxerr_{clr}, bkg_toterr_{clr}, bkg_rmin, bkg_rmax
"""

import numpy as np
import os
from config import *
from my_tools import *


def get_cog(isotab, tflux_max, tfluxerr_max, bkg_r, bkg_toterr):
    """计算单波段的 COG 表。"""
    sma    = isotab['sma'].data
    intens = isotab['intens'].data
    tflux_e = isotab['tflux_e'].data

    # 有效区域：sma < bkg_r 且 intens > 2*bkg_toterr
    ok = (sma < bkg_r) & (intens > 2 * bkg_toterr) & (~np.isnan(intens))
    iso = isotab[ok]

    cog      = iso['tflux_e'] / tflux_max
    cog.name = 'cog'
    cog_lerr = iso['tflux_e'] / (tflux_max + tfluxerr_max)
    cog_lerr.name = 'cog_lerr'
    cog_herr = iso['tflux_e'] / (tflux_max - tfluxerr_max)
    cog_herr.name = 'cog_herr'

    tab = Table()
    tab.add_columns([iso['sma'], iso['intens'], iso['tflux_e'],
                     cog, cog_lerr, cog_herr])
    return tab


def process_chunk(tasks, parms):
    for label, index, clr, row in tasks:
        suffix_clr = f"{clr}_band/{label}_{index}.fits"
        path_in  = os.path.join(parms['iso_dir'], suffix_clr)
        path_out = os.path.join(parms['out_dir'], suffix_clr)

        tflux_max    = row[f'deconv_flux_{clr}']
        tfluxerr_max = row[f'deconv_fluxerr_{clr}']
        bkg_toterr   = row[f'bkg_toterr_{clr}']
        bkg_r = np.sqrt(row['bkg_rmin'] * row['bkg_rmax'])

        isotab = Table.read(path_in)
        tab_cog = get_cog(isotab, tflux_max, tfluxerr_max, bkg_r, bkg_toterr)
        tab_cog.write(path_out)


def run():
    for clr in ['g', 'r', 'i']:
        check_dir(os.path.join(PROPS_ORIGINAL_COG, f"{clr}_band/"), clean=True)
    parms = {'iso_dir': PROPS_ORIGINAL_CLEAN_ISO,
             'out_dir': PROPS_ORIGINAL_COG}
    run_multi_flat(process_chunk, TABLE_PATH, parms)


if __name__ == "__main__":
    run()
