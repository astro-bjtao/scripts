#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成目视检查图。

逻辑：
    将生成的掩模叠在 g/r/i 三色图上，输出 PNG 格式的检查图，
    供后续人工复核使用。
"""

from astropy.table import Table
from astropy.io import fits
import numpy as np
from multiprocessing import Process
import time
import os

# matplotlib 设置（需要在任何画图导入之前设置好 rcParams）
import matplotlib
matplotlib.use('Agg')  # 无 GUI 后端，适合服务器批量出图
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# 路径配置 (所有路径定义在 config.py 里)
from config import *
# 通用工具 (check_dir, get_central, Lognorm2, Band2RGB 等)
from my_tools import *

# ----- 全局 matplotlib 样式设置 -----
plt.rcParams['font.sans-serif'] = "Times New Roman"
plt.rcParams['font.serif'] = "Times New Roman"
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 25
plt.rcParams['axes.labelsize'] = 25
plt.rcParams['figure.figsize'] = (8, 8)


def run_single(table, parms):
    dir_mask_in = parms['mask_in']
    dir_img     = parms['img_dir']
    dir_eyeball = parms['eyeball_dir']

    for i in range(len(table)):
        label = str(table['survey'].data[i], encoding='utf-8')
        index = table['index'].data[i]

        # ---------- 生成掩模 ----------
        suffix = f"{label}_{index}.fits"
        file_mask_in = os.path.join(dir_mask_in, suffix)

        mask_in = fits.getdata(file_mask_in) > 0

        # ---------- 画目视检查图 ----------
        try:
            img_g = fits.getdata(os.path.join(dir_img, f"g_band/{label}_{index}.fits"))
            img_r = fits.getdata(os.path.join(dir_img, f"r_band/{label}_{index}.fits"))
            img_i = fits.getdata(os.path.join(dir_img, f"i_band/{label}_{index}.fits"))

            out_png = os.path.join(dir_eyeball, f"{label}_{index}.png")
            label_text = f"{label.upper()}: {index}"

            plot_eyeball(
                img_g, img_r, img_i,
                out_path=out_png,
                mask=None,
                mode='rgb',          # 三色图 + 蓝掩模
                norm_value=5,
                sw_frac=1,              # 只看中心1/3区域，方便检查
                label_text=label_text
            )
        except FileNotFoundError as e:
            print(f"[WARNING] Missing image for {label}_{index}: {e}")
        except Exception as e:
            print(f"[ERROR] Eyeball plot failed for {label}_{index}: {e}")


def run_multi(path_table, parms, ncpu=120):
    """
    多进程调度：拆表 → 分配进程 → 合并。
    """
    tab = Table.read(path_table)

    n_total = len(tab)
    indices = np.linspace(0, n_total, ncpu + 1, dtype=int)
    process_list = []
    for j in range(ncpu):
        start = indices[j]
        end = indices[j + 1]
        if start == end:
            continue
        sub_tab = tab[start:end]
        process_list.append(
            Process(target=run_single, args=(sub_tab, parms))
        )

    for p in process_list:
        p.start()
        time.sleep(0.01)

    for p in process_list:
        p.join()


def run_all():
    """
    主入口：设置路径，准备输出目录，启动多线程处理。
    """
    check_dir(EYEBALL_IMG_DIR_2, clean=True)

    parms = {
        'mask_in':      TOTAL_MASK_TOTAL,
        'img_dir':      IMG_DIR_2,
        'eyeball_dir':  EYEBALL_IMG_DIR_2      
    }

    run_multi(TABLE_PATH, parms, ncpu=120)


if __name__ == "__main__":
    run_all()