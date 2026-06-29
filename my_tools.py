#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
my_tools.py — 通用天文图像处理与可视化工具集
"""

import os
import math
import shutil
import time
import numpy as np
from multiprocessing import Process
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions
from astropy import units as u
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from photutils.isophote import Isophote, IsophoteList
from photutils.isophote import EllipseGeometry
from photutils.isophote import EllipseSample
from photutils.isophote import Ellipse
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve2d
import matplotlib
matplotlib.use('Agg')  # 无 GUI 后端，适合服务器批量出图
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ============================================================
#  目录与文件操作
# ============================================================
def check_dir(path, clean=False):
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        if clean:
            shutil.rmtree(path)
            os.makedirs(path)

def run_multi(run_single, path_table, parms, ncpu=120):
    """
    通用多进程调度器：拆分星表 → 分配进程 → 等待合并。

    Parameters
    ----------
    run_single : callable
        单进程处理函数，签名为 run_single(table, parms)。
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

# ============================================================
#  图像几何操作
# ============================================================
def get_central(shape):
    ny, nx = shape[0], shape[1]
    cenx = int((nx - 1) / 2)
    ceny = int((ny - 1) / 2)
    return cenx, ceny

def cut_image(data, cenx, ceny, sw):
    """裁切以 (cenx,ceny) 为中心、半径为 sw 的正方形区域（不等尺寸，无 binning）"""
    return data[ceny-sw:ceny+sw+1, cenx-sw:cenx+sw+1]

def bin_image_2x2(data):
    """2x2 像素合并（均值），尺寸减半。用于提升低表面亮度信噪比。"""
    ny, nx = data.shape
    # 裁剪到偶数尺寸
    data = data[:ny//2*2, :nx//2*2]
    return (data[0::2, 0::2] + data[1::2, 0::2] + data[0::2, 1::2] + data[1::2, 1::2]) / 4.0

def cut_bin_image(data, cenx, ceny, sw):
    data = cut_image(data, cenx, ceny, sw)
    return bin_image_2x2(data)

def smooth_img(data, std=1.0):
    """高斯平滑"""
    kernel = Gaussian2DKernel(std).array
    return convolve2d(data, kernel, mode='same', boundary='symm')

# ============================================================
#  形态学操作
# ============================================================
def get_circle(shape, cenx, ceny):

    temp = np.indices(shape, dtype=np.float32)
    xtp = temp[1] - cenx
    ytp = temp[0] - ceny
    im = np.hypot(xtp, ytp)
    return im

def top_hat_kernel(radius):
    """生成半径为 radius 的圆形 top-hat 卷积核"""
    size = 2 * radius + 1
    kernel = get_circle((size, size), radius, radius) < (radius + 0.1)
    return kernel.astype(np.float32)

def extend_mask(mask, ext_r=3):
    """
    对二值掩模进行膨胀操作（卷积 + 阈值）。

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
    生成椭圆形距离矩阵，用于创建椭圆掩模。
    Parameters
    ----------
    shape : tuple (ny, nx)
        图像形状。
    cenx, ceny : int
        椭圆中心坐标。
    b2a : float
        半短轴与半长轴之比 (0, 1]。
    posang : float
        位置角（度），从正X轴逆时针测量。
    Returns
    -------
    im : 2D ndarray
        像素到椭圆中心的归一化距离。等于1时位于椭圆边界上。
    """
    ratio = 1.0 / b2a
    deg2rad = math.pi / 180.0
    cosang = math.cos(deg2rad * posang)
    sinang = math.sin(deg2rad * posang)
    temp = np.indices(shape, dtype=np.float32)
    xtp = temp[1] - cenx
    ytp = temp[0] - ceny
    xtemp = xtp * cosang + ytp * sinang
    ytemp = -xtp * sinang + ytp * cosang
    im = np.hypot(xtemp * ratio, ytemp)
    return im

def reg_to_mask(reg_file, reference_fits, out_mask_fits):
    """
    把 DS9 区域文件转换为二值掩模 FITS 文件。
    
    Parameters
    ----------
    reg_file : str
        输入的 .reg 文件路径。
    reference_fits : str
        参考 FITS 图像路径（用于获取形状和 WCS）。
    out_mask_fits : str
        输出掩模 FITS 文件路径。
    """
    # 读取参考图像获取形状和头部信息
    with fits.open(reference_fits) as hdul:
        ref_data = hdul[0].data
        header = hdul[0].header
    
    shape = ref_data.shape
    
    # 读取区域文件
    regions = Regions.read(reg_file, format='ds9')
    
    # 创建空白掩模 (0 = 背景)
    mask = np.zeros(shape, dtype=np.uint8)
    
    # 将区域转换为像素掩模（需要 WCS 或直接用像素坐标）
    # 如果区域文件使用的是物理坐标，需要 WCS 转换；如果是 image 坐标，则不用。
    # 你给出的 region 文件里写有 "image"，所以是像素坐标，直接转换即可。
    for reg in regions:
        tmp_mask = reg.to_mask(mode='center')
        tmp_mask = tmp_mask.to_image(shape)
        mask[tmp_mask>0]=1
    
    # 保存为 FITS
    fits.PrimaryHDU(data=mask, header=header).writeto(out_mask_fits, overwrite=True)

# ============================================================
#  图像显示与拉伸
# ============================================================
def Lognorm2(data, vmax, a=1000):
    data = data / vmax
    data[data < 1.e-10] = 1.e-10
    data[data > 1.] = 1.
    y = np.log(a * data + 1) / np.log(a + 1)
    y[y > 1.0] = 1.0
    y[y < 0.0] = 0.0
    return y

def Band2RGB(img_r, img_g, img_b):
    assert img_r.shape == img_g.shape == img_b.shape, "Image shapes must be identical"
    img = np.array([img_r, img_g, img_b])
    img = img.transpose((1, 2, 0))
    img[img < 0] = 0
    img[img > 1] = 1
    return img

def plot_scalebar(ax, tot_pix, redshift):
    """
    在图像的 Axes 对象上绘制一个比例尺。

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        目标坐标轴。
    tot_pix : int
        图像的总宽度（像素）。
    pixels_per_kpc : float
        每千秒差距 (kpc) 对应的像素数。
    """
    # 1. 计算图像总宽度的物理尺寸
    kpc_per_pix = kpc_per_pixels(redshift, pix_scale=0.168)
    totwd_kpc = tot_pix * kpc_per_pix

    # 2. 选择一个“整”的比例尺长度 (2, 5, 10, 20...)
    # 这里改用简洁的循环判断
    candidates = [2, 5, 10, 20, 30, 50, 100]
    length_bar = candidates[-1]  # 默认最大值
    for val in candidates:
        if totwd_kpc < val * 4:  # 比例尺长度不超过图像宽度的 1/4
            length_bar = val
            break

    # 3. 计算比例尺在图中的长度
    frac_bar = length_bar / totwd_kpc

    # 4. 在左下角绘制比例尺
    # 横线
    ax.hlines(y=0.1, xmin=0.1, xmax=0.1 + frac_bar, 
              transform=ax.transAxes, color='red', lw=2)
    # 左竖线
    ax.vlines(x=0.1, ymin=0.09, ymax=0.11, 
              transform=ax.transAxes, color='red', lw=2)
    # 右竖线
    ax.vlines(x=0.1 + frac_bar, ymin=0.09, ymax=0.11, 
              transform=ax.transAxes, color='red', lw=2)
    # 文字标注
    ax.text(0.1 + frac_bar/2, 0.13, f'{length_bar} kpc', 
            transform=ax.transAxes, color='red', ha='center', fontsize=12)

# ============================================================
#  通用目视检查图绘制
# ============================================================
def plot_eyeball(
    img_g, img_r, img_i,        # 原始图像 (2D arrays), 可以是 FITS 读取的原始数据
    out_path,                   # 输出 PNG 路径
    mask=None,                  # 可选掩模 (bool 或 01)，将叠加在蓝通道
    mode='rgb',                 # 'rgb', 'rgb_mask', 'deep', 'deep_mask'
    norm_value=3.0,             # 对数拉伸的 vmax
    sw_frac=1.0,                # 裁切范围相对于图像半宽的比例 (0~1)，1=全图, 0.33=中心1/3
    smooth_sigma=1.0,           # 平滑的高斯σ（像素），仅 deep/deep_mask 模式默认平滑
    dpi=300,
    label_text=None              # 显示在左上角的文字，如 "VAGC: 526675"
):
    """
    绘制用于目视检查的星系图像，可选择叠加掩模。
    所有图像均为 2D numpy 数组。
    """
    # 获取中心
    cenx, ceny = get_central(img_g.shape)
    # 计算裁切半宽度
    full_sw = int(min(cenx, ceny))
    sw = int(full_sw * sw_frac)

    # 裁切三个波段图像
    g = cut_image(img_g, cenx, ceny, sw)
    r = cut_image(img_r, cenx, ceny, sw)
    i = cut_image(img_i, cenx, ceny, sw)
    if mask is not None:
        m = cut_image(mask.astype(np.float32), cenx, ceny, sw)
    else:
        m = None

    # 根据模式处理
    if mode in ['deep', 'deep_mask']:
        # 深场模式：先 bin 再平滑
        g = bin_image_2x2(g)
        r = bin_image_2x2(r)
        i = bin_image_2x2(i)
        if m is not None:
            m = bin_image_2x2(m) > 0.5  # 二值化
        # 平滑
        g = smooth_img(g, std=smooth_sigma)
        r = smooth_img(r, std=smooth_sigma)
        i = smooth_img(i, std=smooth_sigma)
        # 取绝对值（处理可能负值）
        g = np.abs(g)
        r = np.abs(r)
        i = np.abs(i)
        # 构建加权灰度图，用于反色背景
        gray = (g + r + i) / 3.0
        # 对数拉伸到 [0,1]
        gray = Lognorm2(gray, vmax=norm_value, a=1000)

        if mode == 'deep':
            # 纯灰度图
            rgb = np.stack([gray, gray, gray], axis=2)
        elif mode == 'deep_mask':
            # 反色背景 + 蓝色掩模
            bg = (1 - gray) * 0.8   # 深色背景
            r_ch = bg
            g_ch = bg
            b_ch = bg
            if m is not None:
                b_ch = b_ch + m * 0.2
            rgb = Band2RGB(r_ch, g_ch, b_ch)
    else:
        # rgb / rgb_mask 模式：不 bin，不平滑，直接对数拉伸
        g = Lognorm2(g, vmax=norm_value, a=1000)
        r = Lognorm2(r, vmax=norm_value, a=1000)
        i = Lognorm2(i, vmax=norm_value, a=1000)
        if mode == 'rgb':
            rgb = Band2RGB(i, r, g)
        elif mode == 'rgb_mask':
            # 减暗一点背景，掩模用蓝色
            g = g * 0.7
            r = r * 0.7
            i = i * 0.7
            if m is not None:
                i = i + m * 0.3
            rgb = Band2RGB(i, r, g)

    # --- 绘制 ---
    fig = plt.figure(figsize=(7, 7), dpi=dpi)
    fig.patch.set_visible(False)
    gs = gridspec.GridSpec(1, 1, figure=fig)
    ax = fig.add_subplot(gs[0])
    ax.imshow(rgb, origin='lower', interpolation='none')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')
    # 移除所有白边
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    plt.margins(0, 0)

    if label_text:
        bbox = {"facecolor": "white", "alpha": 0.9}
        ax.text(0.03, 0.93, label_text, transform=ax.transAxes,
                ha='left', va='top', fontsize=15,
                bbox=bbox, color='black')

    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)


# ============================================================
#  宇宙学工具
# ============================================================
COSMO = FlatLambdaCDM(H0=67.3 * u.km/(u.s*u.Mpc), Tcmb0=2.725, Om0=0.315)

def kpc_per_arcsec(redshift):
    # 使用角直径距离
    distance = COSMO.angular_diameter_distance(redshift).value  # Mpc
    return distance * 4.848 / 1000  # kpc / arcsec

def kpc_per_pixels(redshift, pix_scale=0.168):
    # pix_scale: arcsec / pixel
    kpc_per_arc = kpc_per_arcsec(redshift)
    return pix_scale * kpc_per_arc # kpc / pix

def pixels_per_kpc(redshift, pix_scale=0.168):
    kpc_per_arc = kpc_per_arcsec(redshift)
    return 1.0 / (pix_scale * kpc_per_arc) # pix / kpc

def convert_arcsec2kpc(redshift):
    return kpc_per_arcsec(redshift)  # kpc / arcsec

def Pixs_1kpc(redshift):
    return pixels_per_kpc(redshift, pix_scale=0.168) # pix / kpc

# ============================================================
#  饱和检查专用画图函数
# ============================================================
def plot_eyeball_saturate(
    img_g, img_r, img_i,
    out_path,
    label_text=None,
    sw_frac=1.0,
    dpi=300
):
    """
    生成用于检查星系中心饱和度的三联图：RGB 彩色图 + g-r 颜色图 + g-i 颜色图。
    所有图像均为 2D numpy 数组（原始数据，未拉伸）。
    输出图像无任何白边，三栏无缝拼接。
    """
    # 获取中心
    cenx, ceny = get_central(img_g.shape)

    # 计算裁切半宽度
    full_sw = int(min(cenx, ceny))
    sw = int(full_sw * sw_frac)

    # 裁切图像
    g = cut_image(img_g, cenx, ceny, sw)
    r = cut_image(img_r, cenx, ceny, sw)
    i = cut_image(img_i, cenx, ceny, sw)

    # 计算颜色图（拉伸之前）
    eps = 1.e-10
    color_gr = -2.5 * np.log10(np.maximum(g, eps) / np.maximum(r, eps))
    color_gi = -2.5 * np.log10(np.maximum(g, eps) / np.maximum(i, eps))

    # 中心像素流量
    cen_sw = sw
    cenint_g = g[cen_sw, cen_sw]
    cenint_r = r[cen_sw, cen_sw]
    cenint_i = i[cen_sw, cen_sw]

    # 对数拉伸
    norm_value = np.max([cenint_g, cenint_r, cenint_i])
    g_disp = Lognorm2(g, norm_value, a=1000)
    r_disp = Lognorm2(r, norm_value, a=1000)
    i_disp = Lognorm2(i, norm_value, a=1000)
    rgb_image = Band2RGB(i_disp, r_disp, g_disp)

    # ---- 创建完全无边框的三联图 ----
    fig = plt.figure(figsize=(21, 7), dpi=dpi)
    fig.patch.set_visible(False)

    # 三栏的几何位置：[left, bottom, width, height]，全部填满
    width_per_panel = 1.0 / 3.0
    ax0 = fig.add_axes([0.0, 0.0, width_per_panel, 1.0])
    ax1 = fig.add_axes([width_per_panel, 0.0, width_per_panel, 1.0])
    ax2 = fig.add_axes([2 * width_per_panel, 0.0, width_per_panel, 1.0])

    # 文字样式
    bbox = {"facecolor": "white", "alpha": 1.0}
    styles = {"size": 15, "color": "black", "bbox": bbox}

    # --- 第一栏：RGB 图 ---
    ax0.imshow(rgb_image, origin='lower', interpolation='none', aspect='auto')
    ax0.set_xticks([])
    ax0.set_yticks([])
    ax0.axis('off')
    if label_text:
        ax0.text(0.03, 0.93, label_text, transform=ax0.transAxes,
                 ha='left', linespacing=1.1, **styles)

    # --- 第二栏：g - r 颜色图 ---
    ax1.imshow(color_gr, cmap='cool', vmax=1.5, vmin=-0.5,
               origin='lower', interpolation='none', aspect='auto')
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.axis('off')
    ax1.text(0.05, 0.93, "g - r color image", transform=ax1.transAxes,
             ha='left', linespacing=1.1, **styles)

    # --- 第三栏：g - i 颜色图 + 中心流量 ---
    ax2.imshow(color_gi, cmap='cool', vmax=2, vmin=0,
               origin='lower', interpolation='none', aspect='auto')
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.axis('off')
    ax2.text(0.05, 0.93, "g - i color image", transform=ax2.transAxes,
             ha='left', linespacing=1.1, **styles)

    t = f"central intensity \ng band: {cenint_g:.1f} \nr band: {cenint_r:.1f} \ni band: {cenint_i:.1f}"
    ax2.text(0.05, 0.05, t, transform=ax2.transAxes,
             ha='left', linespacing=1.1, **styles)

    # 保存：不使用 bbox_inches='tight'，直接用固定尺寸，pad_inches=0 确保无白边
    fig.savefig(out_path, dpi=dpi, pad_inches=0.0, bbox_inches=None)
    plt.close(fig)


# ============================================================
#  测量背景
# ============================================================
def Intens2SB(intens, A=27, pixscale=0.168):
    
    return -2.5*np.log10(intens) + A + 5*np.log10(pixscale)


def MeasureBkg(isotab, sma_min):

    # bkg rmin, bkg rmax
    sma_max = sma_min * 1.1**6
    sma     = isotab['sma'].data
    semi_step = 1.1**0.5
    flag = (sma_min/semi_step < sma) & (sma < sma_max*semi_step)
    isotab = isotab[flag]
    # intens, interr, rms of quadrants
    list_intens = []
    list_interr = []
    list_rms = []
    for qx in range(4):
        for i in range(6): # 24 quadrants
            list_intens.append(isotab[f'intens_{qx}'].data[i])
            list_interr.append(isotab[f'intens_err_{qx}'].data[i])
            list_rms.append(isotab[f'rms_{qx}'].data[i])
    # rms and high frequency noise
    rms = np.nanmedian(list_rms)
    err_high = np.nanmedian(list_interr)
    var_high = err_high**2
    # low frequency noise
    std_sky = np.nanstd(list_intens)
    var_low = np.max([0, (std_sky**2 - var_high)])
    var_tot = var_high + var_low
    err_low = np.sqrt(var_low)
    err_tot = np.sqrt(var_tot)
    # local background
    local_bkg = np.nanmedian(list_intens)

    return rms, err_high, err_low, err_tot, local_bkg


# ============================================================
#  Ellipse isophotes
# ============================================================
def convert_wcs(ra, dec, hdr):
    # 0 based
    sky = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    w = WCS(header=hdr)
    xc, yc = w.world_to_pixel(sky)
    xc = float(xc)
    yc = float(yc)

    return xc, yc

def make_masked_image(data, mask):
    # data: science image data
    # mask: source mask data, 1 or True indicate a source pixel need to mask.
    return np.ma.array(data, mask=mask)

def Ellipse_free(masked_image, 
                 x0, y0, eps, pa, 
                 initsma, minsma, maxsma, 
                 step=0.1, fix_center=False, fix_pa=False):
    # masked_image: image data in which soure pixels are masked
    # 0-based x0, y0.
    # eps = 1-b2a
    # pa: anticlockwise from x-axis, unit: degree* np.pi / 180.0
    # fix_center, fix_pa, fix_eps can not be True at same time
    geometry = EllipseGeometry(x0, y0, initsma, eps, pa, 
                               astep=step, linear_growth=False)
    ellipse = Ellipse(masked_image, geometry)
    iso = ellipse.fit_image(sma0=None, 
                            minsma=minsma, 
                            maxsma=maxsma, 
                            step=step,       # default 0.1
                            conver=0.05,     # default 0.05
                            minit=20,
                            maxit=200, 
                            fflag=0.3,
                            maxgerr=0.5,    # default 0.5
                            sclip=3.0, 
                            nclip=0,
                            integrmode='median', 
                            linear=False, 
                            maxrit=None,
                            fix_center=fix_center, 
                            fix_pa=fix_pa, 
                            fix_eps=False)
    
    return iso

def In_Ellipse(masked_image, iso_table):

    isophote_list = []
    len_tab = len(iso_table)
    for i in range(len_tab):
        x0 = iso_table['x0'].data[i]
        y0 = iso_table['y0'].data[i]
        sma = iso_table['sma'].data[i]
        eps = iso_table['ellipticity'].data[i]
        pa = iso_table['pa'].data[i] * np.pi / 180.0
        geometry = EllipseGeometry(x0, y0, sma, eps, pa, 
                                   astep=0.1, linear_growth=False, 
                                   fix_center=True, fix_pa=True, fix_eps=True)
        sample = EllipseSample(masked_image, sma, 
                               sclip=5.0, nclip=0, 
                               integrmode='median',
                               geometry=geometry) # sma need update!
        sample.update(geometry.fix)
        isophote = Isophote(sample, 0, True, stop_code=4)
        isophote_list.append(isophote)
        del x0, y0, sma, eps, pa, geometry

    iso = IsophoteList(isophote_list)

    return iso

def iso_to_table(isophoteList_object):
    # turn to table
    # pa, pa_err turn to degree unit

    return isophoteList_object.to_table(columns='all')

def reshape_isotable(iso_table):
    # iso_table: isophotes in table format
    # reshape_table: reshaped table

    col_names = ['sma', 'intens', 'intens_err', 'ellipticity', 'ellipticity_err', 'pa', 'pa_err', 
                 'x0', 'x0_err', 'y0', 'y0_err', 'rms', 'pix_stddev', 'grad', 'grad_error', 'grad_rerror', 'sarea', 
                 'ndata', 'nflag', 'niter', 'valid', 'stop_code', 'tflux_e', 'tflux_c', 'npix_e', 'npix_c', 
                 'a3', 'b3', 'a4', 'b4', 'a3_err', 'b3_err', 'a4_err', 'b4_err']
    tab = Table()
    for n in col_names:
        dat = iso_table[n].data
        dtp = iso_table[n].dtype
        if dtp == 'object':
            dat = dat.astype(np.float64)
            dtp = np.float64
        col = Column(data=dat, name=n, dtype=dtp)
        tab.add_column(col)
    
    return tab


# ============================================================
#  拟合工具：initsma 查找 + try_fit
# ============================================================

# 等距几何序列，保证不同星系的 sma 对齐
_INIT_SMA_LIST = (300 / 0.168) * 1.1 ** (-1 * np.arange(60))


def get_initsma(target, list_init=None):
    """从几何序列中找最接近 target 的值"""
    if list_init is None:
        list_init = _INIT_SMA_LIST
    delta = np.abs(list_init - target)
    flag = delta < (np.min(delta) + 0.01)
    return float(list_init[flag][0])


def try_fit(masked_image, xc, yc, eps, pa, initsma, minsma, maxsma, step=0.1):
    """尝试自由拟合，成功返回 iso_table，失败返回 None"""
    try:
        iso = Ellipse_free(
            masked_image,
            xc, yc, eps, pa,
            initsma, minsma, maxsma,
            step=step, fix_center=True, fix_pa=False,
        )
        iso_tab = iso_to_table(iso)
        iso_table = reshape_isotable(iso_tab)
        if len(iso_table) > 0:
            return iso_table
    except Exception:
        pass
    return None