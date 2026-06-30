# 处理管线

按顺序执行，每组内可并行。

---

## 0. 基础设施

```
scripts/
  config.py          — 路径配置（所有脚本依赖）
  my_tools/          — 通用工具包
    _io.py           — check_dir, run_multi
    _image.py        — 图像处理
    _plot.py         — 画图（比例尺、目视图）
    _cosmo.py        — 宇宙学换算
    _bkg.py          — 背景测量
    _fitting.py      — Ellipse 拟合 + 矩估计 + 平滑
```

## 1. CUTOUT_SEGMAP — 分割图与目标掩模

| 步骤 | 脚本 | 说明 |
|------|------|------|
| 1.1 | `SE_segmap.py` | SExtractor 分割图 |
| 1.2 | `mask_target.py` | 目标星系掩模 |
| 1.3 | `eyeball_mask.py` | 目视检查掩模 |
| 1.4 | `eyeball_mask_target.py` | 目视检查目标掩模 |

## 2. MASK — 合成总掩模 `mask/`

| 步骤 | 脚本 | 说明 |
|------|------|------|
| 2.1 | `mask/mask_outer.py` | 外区掩模（椭圆环） |
| 2.2 | 手动编辑 | DS9 微调 |
| 2.3 | `mask/mask_inner.py` | 内区掩模 |
| 2.4 | `mask/mask_segmap_all.py` | 合并分割图掩模 |
| 2.5 | `mask/mask_star.py` | 亮星掩模 |
| 2.6 | `mask/mask_companion.py` | 伴星系掩模 |
| 2.7 | `mask/mask_total_auto.py` | 自动合成 |
| 2.8 | `mask/mask_manual.py` | 手动修正 |
| 2.9 | `mask/mask_total.py` | 最终总掩模 → `TOTAL_MASK_TOTAL` |

## 3. BACKGROUND — 背景测量与图像预处理 `bkg/`

| 步骤 | 脚本 | 说明 |
|------|------|------|
| 3.1 | `bkg/profile_quarters.py` | 四象限圆形孔径测光 → `LIMIT_DEPTH_ISOTAB` |
| 3.2 | `bkg/get_bkg_rmin.py` | 确定背景区域起止半径 |
| 3.3 | `bkg/measure_bkg_err.py` | 测量背景误差、极限面亮度 |
| 3.4 | `bkg/subtract_bkg.py` | 减去局部背景 → `IMG_DIR_2` |
| 3.5 | `bkg/sum_all.py` | g+r+i 合成 a 波段 |
| 3.6 | `bkg/eyeball_bkg.py` | 目视检查背景区域 |

## 4. FITTING — 等强度线拟合 `fitting/`

| 步骤 | 脚本 | 说明 |
|------|------|------|
| 4.1 | `fitting/free_fitting_original.py` | Ellipse 拟合 + moments_estimate + smooth_isotable → `PROPS_ORIGINAL_ISOTAB` |
| 4.2 | `fitting/eyeball_fitting.py` | 1Re/2Re 目视检查图 → `PROPS_ORIGINAL_EYEBALL_FITTING` |

## 5. EYEBALL — 辅助目视检查 `eyeball/`

| 步骤 | 脚本 | 说明 |
|------|------|------|
| 5.1 | `eyeball/eyeball.py` | 掩模叠加目视检查 |
| 5.2 | `eyeball/eyeball_saturate.py` | 饱和检查三联图 |

---

运行示例：
```bash
cd /data1/bjtao/StellarHalo_z02/Process2/scripts
python mask/mask_total.py
python bkg/profile_quarters.py
python fitting/free_fitting_original.py
python fitting/eyeball_fitting.py
```
