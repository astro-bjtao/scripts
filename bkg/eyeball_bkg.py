#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
画出背景区域，并列出极限星等。
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


### Run it

def run_single(table, parms):

    # parms
    dir_img     = parms['img_dir']
    dir_eyeball = parms['eyeball_dir']

    for i in range(len(table)):
        # information
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]
        bkg_rmin = table['bkg_rmin'].data[i]
        bkg_rmax = table['bkg_rmax'].data[i]
        limit_mu_g = table[f'limit_mu_g'].data[i]
        limit_mu_r = table[f'limit_mu_r'].data[i]
        limit_mu_i = table[f'limit_mu_i'].data[i]
        # path
        path_img_g = os.path.join(dir_img, f"g_band/{label}_{index}.fits")
        path_img_r = os.path.join(dir_img, f"r_band/{label}_{index}.fits")
        path_img_i = os.path.join(dir_img, f"i_band/{label}_{index}.fits")
        path_eyeball = os.path.join(dir_eyeball, f"{label}_{index}.png")
        # read, image
        image_g = fits.getdata(path_img_g)
        image_r = fits.getdata(path_img_r)
        image_i = fits.getdata(path_img_i)
        # cut and bin image
        cenx, ceny = get_central(image_g.shape)
        #norm_value = np.max([image_g[ceny,cenx], image_r[ceny,cenx], image_i[ceny,cenx]])/2
        norm_value = 3
        sw = int(np.min([cenx, ceny]))
        # cut image
        image_g = cut_image(image_g, cenx, ceny, sw)
        image_r = cut_image(image_r, cenx, ceny, sw)
        image_i = cut_image(image_i, cenx, ceny, sw)
        # LogStretch
        image_g = Lognorm2(image_g, norm_value, a=1000)
        image_r = Lognorm2(image_r, norm_value, a=1000)
        image_i = Lognorm2(image_i, norm_value, a=1000)
        # rgb image
        rgb_image = Band2RGB(image_i, image_r, image_g)

        # Create drawing objects
        fig = plt.figure(figsize=(7,7),dpi=300) # constrained_layout=True
        plt.subplots_adjust(top=1, bottom=0,right=1,left=0, hspace=0, wspace=0)
        plt.margins(0, 0)
        gs = gridspec.GridSpec(1, 1, figure=fig)

        plt.xticks([])
        plt.yticks([])

        # Create main image area
        ax0 = fig.add_subplot(gs[0])
        
        # show image
        ax0.imshow(rgb_image, origin='lower', interpolation='none')
        ax0.set_xticks([])
        ax0.set_yticks([])
        ax0.axis('off')

        # bkg region
        bkg_rmin /= 1.1**0.5
        bkg_rmax *= 1.1**0.5
        for sma in [bkg_rmin, bkg_rmax]:
            ellipse = patches.Ellipse((sw, sw), 2*sma,  2*sma, angle=0, alpha=1.0, 
                                    fill=False, color='red', linewidth=2, linestyle='--')
            ax0.add_artist(ellipse)
        
        # information text
        bbox = {"facecolor": "white", "alpha": 1.0}
        styles = {"size": 25, "color": "black", "bbox": bbox}
        t = f"{label.upper()}: {index}"
        ax0.text(0.03,0.93, t, transform=ax0.transAxes, ha='left', linespacing=1.1, **styles)

        t = f"limit_mu_g: {limit_mu_g:.3f} \nlimit_mu_r: {limit_mu_r:.3f} \nlimit_mu_i: {limit_mu_i:.3f}"
        ax0.text(0.03,0.03, t, transform=ax0.transAxes, ha='left', linespacing=1.1, **styles)
        
        # save figure
        fig.savefig(path_eyeball,  bbox_inches="tight", pad_inches=0.1)
        plt.close()

def run_all():
    """
    主入口：设置路径，准备输出目录，启动多线程处理。
    """
    check_dir(LIMIT_DEPTH_EYEBALL_BKG, clean=True)

    parms = {
        'img_dir':      IMG_DIR,
        'eyeball_dir':  LIMIT_DEPTH_EYEBALL_BKG      
    }

    run_multi(run_single, TABLE_PATH, parms, ncpu=120)

if __name__ == "__main__":

    run_all()