# StellarHalo Process2 — 星系测光管线

对 `StellarHalo` 巡天 z~0.02 星系样本进行图像预处理、掩模生成、背景测量和等强度线拟合的自动化脚本集。

## 目录结构

```
scripts/
  config.py             路径配置（所有脚本依赖）
  my_tools/             通用工具包（7 模块）
    _io.py              check_dir, run_multi, run_multi_flat
    _image.py           图像几何/形态/拉伸
    _cosmo.py           宇宙学距离换算
    _bkg.py             背景测量
    _fitting.py         Ellipse 拟合、矩估计、等强度线平滑
    _plot.py            比例尺、目视检查图
  mask/                 掩模生成（8 脚本）
  bkg/                  背景测量 + 图像预处理（6 脚本）
  fitting/              等强度线拟合 + 建模 + 清理（4 脚本）
  eyeball/              辅助目视检查（2 脚本）
  compare_old/          重构前后对比档案
  PIPELINE.md           处理管线顺序文档
  LOG.md                工作日志
```

## 安装

在服务器 `/data1/bjtao/miniconda3/lib/python3.9/site-packages/` 下放置 `proc2.pth`：

```
/data1/bjtao/StellarHalo_z02/Process2/scripts
```

之后任何位置可 `import config` / `import my_tools`。

## 完整管线

```bash
cd /data1/bjtao/StellarHalo_z02/Process2/scripts

# 1. 掩模
python mask/mask_total.py

# 2. 背景测量
python bkg/profile_quarters.py
python bkg/get_bkg_rmin.py
python bkg/measure_bkg_err.py

# 3. 图像预处理
python bkg/subtract_bkg.py
python bkg/sum_all.py

# 4. 等强度线拟合
python fitting/free_fitting_original.py      # 自动矩估计 + 平滑, 98/98

# 5. 椭圆建模 + 像素清理
python fitting/build_model.py                # 科学图像
python fitting/build_model.py var            # 方差图像
python fitting/clean_pixels.py               # 生成干净图像
python fitting/clean_pixels.py var           # 生成干净方差

# 6. 目视检查
python fitting/eyeball_fitting.py
```

## 核心算法

| 函数 | 模块 | 说明 |
|------|------|------|
| `moments_estimate` | `_fitting` | 通量加权二阶矩自动估计 eps/PA 初值 |
| `smooth_isotable` | `_fitting` | SG 滤波 + sigma-clip + PA 解缠 + 漂移检测，修复跃变/重叠 |
| `try_fit` | `_fitting` | 带异常捕获的封装 |
| `run_multi_flat` | `_io` | 星系×波段展开，跑满多核 |
| `plot_eyeball` | `_plot` | RGB 三色目视检查图 |

## 依赖

```
astropy  numpy  scipy  matplotlib  photutils  regions
```

## 数据

图像和星表路径由 `config.py` 统一管理。数据根目录：`/data1/bjtao/StellarHalo_z02/Process2/`
