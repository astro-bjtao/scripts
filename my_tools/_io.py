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


def run_multi_flat(run_chunk, path_table, parms, bands=None, ncpu=120):
    """
    通用多进程调度器：星系×波段展开后均分到 ncpu 个进程。

    Parameters
    ----------
    run_chunk : callable
        单进程处理函数，签名为 ``run_chunk(tasks, parms)``。
        每个 task 为 (label_str, index_int, clr_str, row_dict)。
    path_table : str
        主星表 FITS 文件路径。
    parms : dict
        传递给 run_chunk 的参数字典。
    bands : list of str
        波段列表，默认 ['g', 'r', 'i']。
    ncpu : int
        并行进程数，默认 120。
    """
    from astropy.table import Table

    if bands is None:
        bands = ['g', 'r', 'i']

    tab = Table.read(path_table)

    # 展开星系×波段，每个 task 包含 row dict 方便读取其他参数
    all_tasks = []
    for row in tab:
        # 兼容 bytes 和 str 类型
        val = row['survey']
        label = val.decode('utf-8') if isinstance(val, bytes) else str(val)
        index = row['index']
        row_dict = {col: row[col] for col in row.colnames}
        for clr in bands:
            all_tasks.append((label, index, clr, row_dict))

    n_total = len(all_tasks)
    indices = np.linspace(0, n_total, ncpu + 1, dtype=int)

    processes = []
    for j in range(ncpu):
        start, end = indices[j], indices[j + 1]
        if start == end:
            continue
        processes.append(
            Process(target=run_chunk, args=(all_tasks[start:end], parms))
        )

    for p in processes:
        p.start()
        time.sleep(0.01)

    for p in processes:
        p.join()
