#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
批量生成星系饱和检查图。

对每个星系，输出一张三联图：
    左：RGB 三色图
    中：g - r 颜色图
    右：g - i 颜色图（同时标注中心像素流量）
"""

from astropy.table import Table
from astropy.io import fits
import numpy as np
import os

# 路径配置
from config import *
# 通用工具
from my_tools import check_dir, plot_eyeball_saturate, run_multi


def run_single(table, parms):
    """
    单线程处理函数：遍历星表，对每个星系生成饱和检查图。
    """
    dir_img = parms['img_dir']
    dir_eyeball = parms['eyeball_dir']

    for i in range(len(table)):
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]

        # 图像路径
        file_g = os.path.join(dir_img, f"g_band/{label}_{index}.fits")
        file_r = os.path.join(dir_img, f"r_band/{label}_{index}.fits")
        file_i = os.path.join(dir_img, f"i_band/{label}_{index}.fits")
        out_file = os.path.join(dir_eyeball, f"{label}_{index}.png")

        # 跳过缺失图像
        if not (os.path.exists(file_g) and os.path.exists(file_r) and os.path.exists(file_i)):
            print(f"[WARNING] Missing image(s) for {label}_{index}, skipping.")
            continue

        # 读取图像
        img_g = fits.getdata(file_g)
        img_r = fits.getdata(file_r)
        img_i = fits.getdata(file_i)

        # 画图
        plot_eyeball_saturate(
            img_g, img_r, img_i,
            out_path=out_file,
            label_text=f"{label.upper()}: {index}",
            sw_frac=0.2,              # 全图裁切；如需只看中心可以改小
            dpi=300
        )


def run_all():
    """
    主入口：设置路径，准备输出目录，启动多进程处理。
    """
    check_dir(EYEBALL_SATURATE_DIR, clean=True)

    parms = {
        'img_dir': IMG_DIR_2,   # 或者你的图像目录，按需修改
        'eyeball_dir': EYEBALL_SATURATE_DIR
    }

    run_multi(run_single, TABLE_PATH, parms, ncpu=64)


if __name__ == "__main__":
    run_all()
