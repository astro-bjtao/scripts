#!/usr/bin/env python3
"""
功能对比审计：原始 my_tools.py vs 重构 my_tools/ 包。

用法：
    cd /data1/bjtao/StellarHalo_z02/Process2/scripts
    python3 compare_old/audit_comparison.py

结果（2026-06-30）：55 passed, 0 failed — 重构前后完全等价。
"""

import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# 加载原始版本
orig_ns = {}
exec(open(os.path.join(os.path.dirname(__file__), 'my_tools_original.py')).read(), orig_ns)

# 加载重构版本
from my_tools._image import *
from my_tools._cosmo import *
from my_tools._bkg import *
from my_tools._fitting import *
from my_tools._fitting import _INIT_SMA_LIST

old = orig_ns
ok = fail = 0

def test(name, orig, new, rtol=1e-10):
    global ok, fail
    try:
        if isinstance(orig, np.ndarray):
            match = np.allclose(orig, new, rtol=rtol)
        elif isinstance(orig, float):
            match = np.isclose(orig, new, rtol=rtol)
        elif orig is None:
            match = new is None
        else:
            match = (orig == new)
        if match: ok += 1
        else: fail += 1; print(f"  MISMATCH {name}")
    except Exception as e:
        fail += 1; print(f"  ERROR {name}: {e}")

# Geometry
d = np.arange(100).reshape(10,10).astype(float)
test("get_central",      old['get_central']((10,20)), get_central((10,20)))
test("cut_image",        old['cut_image'](d,5,5,3), cut_image(d,5,5,3))
test("bin_image_2x2",    old['bin_image_2x2'](d), bin_image_2x2(d))
test("cut_bin_image",    old['cut_bin_image'](d,5,5,3), cut_bin_image(d,5,5,3))
d2 = np.random.RandomState(42).randn(20,20)
test("smooth_img",       old['smooth_img'](d2,1.0), smooth_img(d2,1.0))
test("get_circle",       old['get_circle']((10,10),5,5), get_circle((10,10),5,5))
test("top_hat_kernel",   old['top_hat_kernel'](3), top_hat_kernel(3))
m = np.zeros((20,20),bool); m[10,10]=True
test("extend_mask",      old['extend_mask'](m,2), extend_mask(m,2))
test("get_ellipse",      old['get_ellipse']((10,10),5,5,0.5,45), get_ellipse((10,10),5,5,0.5,45))

# Display
d3 = np.array([0.01,0.5,1.5,2.0])
test("Lognorm2",         old['Lognorm2'](d3.copy(),1.0), Lognorm2(d3.copy(),1.0))
a,b,c = np.eye(3),2*np.eye(3),0.5*np.eye(3)
test("Band2RGB",         old['Band2RGB'](a,b,c), Band2RGB(a,b,c))
m2 = np.zeros((5,5),bool); m2[2,2]=True
test("make_masked_image",old['make_masked_image'](np.ones((5,5)),m2).mask.sum(), make_masked_image(np.ones((5,5)),m2).mask.sum())

# Cosmo
test("kpc_per_arcsec",     old['kpc_per_arcsec'](0.02), kpc_per_arcsec(0.02))
test("kpc_per_pixels",     old['kpc_per_pixels'](0.02), kpc_per_pixels(0.02))
test("pixels_per_kpc",     old['pixels_per_kpc'](0.02), pixels_per_kpc(0.02))
test("convert_arcsec2kpc", old['convert_arcsec2kpc'](0.02), convert_arcsec2kpc(0.02))
test("Pixs_1kpc",          old['Pixs_1kpc'](0.02), Pixs_1kpc(0.02))

# Bkg
test("Intens2SB", old['Intens2SB'](1.0), Intens2SB(1.0))

# Fitting
test("get_initsma", old['get_initsma'](100.0), get_initsma(100.0))
test("_INIT_SMA_LIST", old['_INIT_SMA_LIST'], _INIT_SMA_LIST)

# reshape_isotable
from astropy.table import Table, Column
cols = ["sma","intens","intens_err","ellipticity","ellipticity_err","pa","pa_err",
        "x0","x0_err","y0","y0_err","rms","pix_stddev","grad","grad_error",
        "grad_rerror","sarea","ndata","nflag","niter","valid","stop_code",
        "tflux_e","tflux_c","npix_e","npix_c","a3","b3","a4","b4",
        "a3_err","b3_err","a4_err","b4_err"]
t = Table()
for cn in cols:
    t.add_column(Column(data=[1.,2.,3.], name=cn, dtype=np.float64))
t['stop_code'] = [0,0,4]
try:
    o_r = old['reshape_isotable'](t)
    n_r = reshape_isotable(t)
    for cn in cols:
        test(f"reshape_isotable[{cn}]", o_r[cn].data, n_r[cn].data)
except Exception as e:
    fail += 1; print(f"  ERROR reshape_isotable: {e}")

# check_dir
import tempfile, shutil
td = tempfile.mkdtemp()
old['check_dir'](td+"/test_orig")
from my_tools._io import check_dir as new_check_dir
new_check_dir(td+"/test_new")
test("check_dir_exist", os.path.isdir(td+"/test_orig"), os.path.isdir(td+"/test_new"))
shutil.rmtree(td)

print(f"\n{'='*50}")
print(f"  {ok} passed, {fail} failed")
if fail == 0:
    print("  ALL FUNCTIONS VERIFIED IDENTICAL")
