# StellarHalo Process2 — 星系测光管线

对 `StellarHalo` 巡天 z~0.02 星系样本进行图像预处理、掩模生成、背景测量和等强度线拟合的一套自动化脚本。

## 目录结构

```
scripts/
  config.py             路径配置
  my_tools/             通用工具包
    _io.py              目录操作、多进程调度
    _image.py           图像几何变换、形态学、显示拉伸
    _cosmo.py           宇宙学距离换算
    _bkg.py             背景测量
    _fitting.py         Ellipse 等强度线拟合、矩估计、平滑后处理
    _plot.py            比例尺、RGB 目视图、饱和检查三联图
  mask/                 掩模生成（8 脚本）
  bkg/                  背景测量与图像预处理（6 脚本）
  fitting/              等强度线拟合与目视检查（2 脚本）
  eyeball/              辅助目视检查（2 脚本）
  PIPELINE.md           处理管线顺序文档
```

## 安装

在服务器 `/data1/bjtao/miniconda3/lib/python3.9/site-packages/` 下放置 `proc2.pth`：

```
/data1/bjtao/StellarHalo_z02/Process2/scripts
```

之后任何位置可 `import config` / `import my_tools`，Jupyter notebook 也可直接调用。

## 快速开始

```bash
cd /data1/bjtao/StellarHalo_z02/Process2/scripts

# 掩模
python mask/mask_total.py

# 背景测量
python bkg/profile_quarters.py
python bkg/subtract_bkg.py

# 等强度线拟合（自动矩估计 + 平滑）
python fitting/free_fitting_original.py

# 目视检查
python fitting/eyeball_fitting.py
```

## 核心算法

| 函数 | 说明 |
|------|------|
| `moments_estimate` | 通量加权二阶矩自动估计 eps/PA 初值 |
| `smooth_isotable` | SG 滤波 + sigma-clip + 插值，修复等强度线跃变和轮廓重叠 |
| `try_fit` | 带异常捕获的自由 Ellipse 拟合封装 |

## 依赖

```
astropy  numpy  scipy  matplotlib  photutils  regions
```

## 数据

图像和星表路径由 `config.py` 统一管理。当前数据根目录：`/data1/bjtao/StellarHalo_z02/Process2/`
