#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
目录与文件操作、多进程调度。
"""

import os
import shutil
import time
import numpy as np
from multiprocessing import Process


def check_dir(path, clean=False):
    """确保目录存在。若 clean=True 则先清空重建。"""
    if clean and os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


def run_multi(run_single, path_table, parms, ncpu=120):
    """
    通用多进程调度器：拆分星表 → 分配进程 → 等待合并。

    Parameters
    ----------
    run_single : callable
        单进程处理函数，签名为 ``run_single(table, parms)``。
    path_table : str
        主星表 FITS 文件路径。
    parms : dict
        传递给 run_single 的参数字典。
    ncpu : int
        并行进程数，默认 120。
    """
    from astropy.table import Table

    tab = Table.read(path_table)
    n_total = len(tab)
    indices = np.linspace(0, n_total, ncpu + 1, dtype=int)

    processes = []
    for j in range(ncpu):
        start, end = indices[j], indices[j + 1]
        if start == end:
            continue
        processes.append(
            Process(target=run_single, args=(tab[start:end], parms))
        )

    for p in processes:
        p.start()
        time.sleep(0.01)

    for p in processes:
        p.join()
