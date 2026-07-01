#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ellipse 等强度线拟合、矩估计初值、等强度线后处理平滑。
"""

import numpy as np
from astropy import units as u
from astropy.wcs import WCS
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from photutils.isophote import (
    Isophote, IsophoteList,
    EllipseGeometry, EllipseSample, Ellipse,
)

# ============================================================
#  WCS 坐标转换
# ============================================================

def convert_wcs(ra, dec, hdr):
    """
    从 WCS 将 (RA, Dec) 转为像素坐标 (xc, yc)，0-indexed。

    Parameters
    ----------
    ra, dec : float
        赤经、赤纬（度）。
    hdr : astropy.io.fits.Header
        FITS 头，需含 WCS 关键字。

    Returns
    -------
    xc, yc : float
    """
    sky = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    w = WCS(header=hdr)
    xc, yc = w.world_to_pixel(sky)
    return float(xc), float(yc)


# ============================================================
#  Ellipse 自由拟合
# ============================================================

def Ellipse_free(masked_image,
                 x0, y0, eps, pa,
                 initsma, minsma, maxsma,
                 step=0.1, fix_center=False, fix_pa=False):
    """
    对掩模图像执行自由 Ellipse 拟合。

    Parameters
    ----------
    masked_image : np.ma.MaskedArray
    x0, y0 : float
        初始中心（0-indexed）。
    eps : float
        初始椭圆率 1 - b/a。
    pa : float
        初始位置角（弧度），从 X+ 逆时针。
    initsma : float
        起始半长轴（像素）。
    minsma, maxsma : float
        最小/最大半长轴。
    step : float
        径向步长（对数增长因子），默认 0.1。
    fix_center, fix_pa : bool
        是否固定中心/PA。

    Returns
    -------
    iso : photutils.isophote.IsophoteList
    """
    geometry = EllipseGeometry(x0, y0, initsma, eps, pa,
                               astep=step, linear_growth=False)
    ellipse = Ellipse(masked_image, geometry)
    iso = ellipse.fit_image(
        sma0=None, minsma=minsma, maxsma=maxsma,
        step=step, conver=0.05, minit=20, maxit=200,
        fflag=0.3, maxgerr=0.5, sclip=3.0, nclip=0,
        integrmode='median', linear=False, maxrit=None,
        fix_center=fix_center, fix_pa=fix_pa, fix_eps=False,
    )
    return iso


# ============================================================
#  结果格式转换
# ============================================================

def In_Ellipse(masked_image, iso_table):
    """根据 iso_table 重建固定几何的 IsophoteList（用于重采样）。"""
    isophote_list = []
    for row in iso_table:
        x0, y0 = row['x0'], row['y0']
        sma = row['sma']
        eps = row['ellipticity']
        pa = row['pa'] * np.pi / 180.0
        geometry = EllipseGeometry(x0, y0, sma, eps, pa,
                                   astep=0.1, linear_growth=False,
                                   fix_center=True, fix_pa=True, fix_eps=True)
        sample = EllipseSample(masked_image, sma,
                               sclip=5.0, nclip=0, integrmode='median',
                               geometry=geometry)
        sample.update(geometry.fix)
        isophote_list.append(Isophote(sample, 0, True, stop_code=4))
    return IsophoteList(isophote_list)


def iso_to_table(isophote_list):
    """IsophoteList → astropy Table（全部列，PA 转为度）。"""
    return isophote_list.to_table(columns='all')


def reshape_isotable(iso_table):
    """
    标准化等强度线表：提取指定列，object 类型转为 float64。

    确保下游代码访问列时不因类型问题报错。
    """
    col_names = [
        'sma', 'intens', 'intens_err',
        'ellipticity', 'ellipticity_err', 'pa', 'pa_err',
        'x0', 'x0_err', 'y0', 'y0_err',
        'rms', 'pix_stddev', 'grad', 'grad_error', 'grad_rerror',
        'sarea', 'ndata', 'nflag', 'niter', 'valid', 'stop_code',
        'tflux_e', 'tflux_c', 'npix_e', 'npix_c',
        'a3', 'b3', 'a4', 'b4',
        'a3_err', 'b3_err', 'a4_err', 'b4_err',
    ]
    tab = Table()
    for name in col_names:
        col_dtype = iso_table[name].dtype        # 列级别的 dtype
        data = np.array(iso_table[name].data)    # 统一转为 ndarray（处理 memoryview）
        if col_dtype == object or data.dtype == object:
            data = data.astype(np.float64)
            col_dtype = np.float64
        tab.add_column(Column(data=data, name=name, dtype=col_dtype))
    return tab


# ============================================================
#  拟合辅助
# ============================================================

# 对数等距几何序列，保证不同星系 sma 对齐
_INIT_SMA_LIST = (300.0 / 0.168) * 1.1 ** (-np.arange(60))


def get_initsma(target, list_init=None):
    """从几何序列中取最接近 target 的值。"""
    if list_init is None:
        list_init = _INIT_SMA_LIST
    delta = np.abs(list_init - target)
    best = delta < (delta.min() + 0.01)
    return float(list_init[best][0])


def try_fit(masked_image, xc, yc, eps, pa, initsma, minsma, maxsma, step=0.1):
    """
    尝试自由拟合；成功返回 iso_table，失败返回 None。

    使用 fix_center=True, fix_pa=False。
    """
    try:
        iso = Ellipse_free(masked_image, xc, yc, eps, pa,
                           initsma, minsma, maxsma,
                           step=step, fix_center=True, fix_pa=False)
        tab = iso_to_table(iso)
        iso_table = reshape_isotable(tab)
        if len(iso_table) > 0:
            return iso_table
    except Exception:
        pass
    return None


# ============================================================
#  图像矩估计：自动估计 eps / PA 初值
# ============================================================

def moments_estimate(masked_image, xc, yc, r_in, r_out):
    """
    在环形区域内计算通量加权二阶矩，估计椭圆率 eps 和位置角 PA。

    原理
    ----
    通量加权的二阶中心矩构成对称矩阵 Q：

        Q_xx = Σ I (x - x̄)² / Σ I
        Q_yy = Σ I (y - ȳ)² / Σ I
        Q_xy = Σ I (x - x̄)(y - ȳ) / Σ I

    其中 (x̄, ȳ) 是环形区域内的通量质心。对 Q 做特征分解：

        a² = (Q_xx + Q_yy)/2 + √[((Q_xx-Q_yy)/2)² + Q_xy²]
        b² = (Q_xx + Q_yy)/2 - √[((Q_xx-Q_yy)/2)² + Q_xy²]

    eps = 1 - b/a  （photutils 约定，0 = 正圆）
    PA  = ½ atan2(2 Q_xy, Q_xx - Q_yy)  （X+ 逆时针，度）

    Parameters
    ----------
    masked_image : np.ma.MaskedArray
    xc, yc : float
        初始中心坐标（0-indexed）。
    r_in, r_out : float
        环形区域内/外半径（像素）。

    Returns
    -------
    eps : float or None
    pa_deg : float or None
        [0°, 180°)，X+ 逆时针。可直接用于 EllipseGeometry。
    """
    ny, nx = masked_image.data.shape
    y, x = np.indices((ny, nx), dtype=np.float64)

    r = np.hypot(x - xc, y - yc)
    ring = (r > r_in) & (r < r_out) & (~masked_image.mask)
    if ring.sum() < 10:
        return None, None

    flux = masked_image.data[ring]
    xr, yr = x[ring], y[ring]
    total = flux.sum()
    x_c = (flux * xr).sum() / total
    y_c = (flux * yr).sum() / total

    dx, dy = xr - x_c, yr - y_c
    qxx = (flux * dx * dx).sum() / total
    qyy = (flux * dy * dy).sum() / total
    qxy = (flux * dx * dy).sum() / total

    d = np.sqrt(((qxx - qyy) / 2) ** 2 + qxy ** 2)
    a2 = (qxx + qyy) / 2 + d
    b2 = (qxx + qyy) / 2 - d
    if a2 <= 0 or b2 <= 0 or b2 > a2:
        return None, None

    eps = float(1.0 - np.sqrt(b2 / a2))
    pa_rad = 0.5 * np.arctan2(2.0 * qxy, qxx - qyy)
    pa_deg = pa_rad * 180.0 / np.pi
    if pa_deg < 0:
        pa_deg += 180.0

    return eps, float(pa_deg)


# ============================================================
#  等强度线后处理：平滑跃变、修复重叠
# ============================================================

def smooth_isotable(iso, sigma_clip=3.0, window=7):
    """
    平滑等强度线的 eps 和 PA，修复跃变和轮廓重叠。

    算法
    ----
    1. PA 相位解缠：消除 0°/180° 边界虚假跃变
    2. 检测 PA 系统性漂移（范围 > 45°），自适应扩大 SG 窗口
    3. log(sma) 空间 SG 滤波 → sigma-clip → 异常点 log-linear 插值
    4. 确保半短轴单调递增，消除轮廓重叠

    Parameters
    ----------
    iso : Table
        等强度线表（需 sma, ellipticity, pa, stop_code 列）。
    sigma_clip : float
        异常阈值（sigma），默认 3.0。
    window : int
        SG 滤波器窗口（奇数），默认 7。

    Returns
    -------
    iso_fixed : Table
    stats : dict
        n_fixed_eps, n_fixed_pa, n_fixed_overlap
    """
    from scipy.signal import savgol_filter

    iso = iso.copy()
    n = len(iso)
    sma = iso['sma'].data.copy()
    eps = iso['ellipticity'].data.copy()
    pa = iso['pa'].data.copy()
    stop = iso['stop_code'].data.copy()
    log_s = np.log(sma)

    stats = {'n_fixed_eps': 0, 'n_fixed_pa': 0, 'n_fixed_overlap': 0}

    # ---------- PA 相位解缠 ----------
    def _unwrap(pa_arr):
        pu = pa_arr.copy().astype(np.float64)
        for i in range(1, len(pu)):
            d = pu[i] - pu[i - 1]
            if d > 90:
                pu[i:] -= 180
            elif d < -90:
                pu[i:] += 180
        return pu

    def _rewrap(pu):
        return pu % 180.0

    # Step 0：stop_code==4 继承前一条有效值
    for i in range(1, n):
        if stop[i] == 4:
            eps[i] = eps[i - 1]
            pa[i] = pa[i - 1]
            stats['n_fixed_eps'] += 1
            stats['n_fixed_pa'] += 1

    valid = (stop != 4)
    nv = valid.sum()
    if nv < 5:
        iso['ellipticity'] = eps
        iso['pa'] = pa
        return iso, stats

    # 窗口自适应
    w = window | 1  # 确保奇数
    w = max(3, min(w, nv - 1 + (nv % 2)))

    idx = np.where(valid)[0]
    order = np.argsort(log_s[valid])
    xs = log_s[idx][order]
    es = eps[idx][order]
    ps_raw = pa[idx][order]

    ps = _unwrap(ps_raw)

    # PA 漂移检测 → 大窗口
    w_pa = w
    if np.ptp(ps) > 45:
        w_pa = max(w * 2, (nv // 3) | 1)
        w_pa = min(w_pa, nv - 1 + (nv % 2))

    es_sg = savgol_filter(es, w, 2)
    ps_sg = savgol_filter(ps, w_pa, 2)

    # Sigma-clip
    r_eps = np.abs(es - es_sg)
    r_pa = np.abs(ps - ps_sg)
    th_eps = sigma_clip * np.median(r_eps) * 1.4826 + 0.02
    th_pa = sigma_clip * np.median(r_pa) * 1.4826 + 3.0

    bad = (r_eps > th_eps) | (r_pa > th_pa)
    stats['n_fixed_eps'] += (r_eps > th_eps).sum()
    stats['n_fixed_pa'] += (r_pa > th_pa).sum()

    # 插值修复
    if bad.any() and (~bad).sum() >= 2:
        good = ~bad
        es[bad] = np.interp(xs[bad], xs[good], es[good])
        ps[bad] = np.interp(xs[bad], xs[good], ps[good])

    eps_fin = eps.copy()
    pa_fin = pa.copy()
    eps_fin[idx[order]] = es
    pa_fin[idx[order]] = _rewrap(ps)

    # 半短轴单调性
    b = sma * (1.0 - eps_fin)
    for i in range(1, n):
        if b[i] <= b[i - 1]:
            b_min = b[i - 1] * 1.001
            eps_fin[i] = max(0.0, min(0.99, 1.0 - b_min / sma[i]))
            b[i] = sma[i] * (1.0 - eps_fin[i])
            stats['n_fixed_overlap'] += 1

    iso['ellipticity'] = eps_fin
    iso['pa'] = pa_fin
    return iso, stats
