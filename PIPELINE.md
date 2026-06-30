# 处理管线

按顺序执行，每组内可并行。

---

## 1. CUTOUT_SEGMAP — 分割图与目标掩模

| 步骤 | 脚本 | 说明 |
|------|------|------|
| 1.1 | `SE_segmap.py` | SExtractor 生成分割图 |
| 1.2 | `mask_target.py` | 目标星系掩模 |
| 1.3 | `eyeball_mask.py` | 目视检查掩模 |
| 1.4 | `eyeball_mask_target.py` | 目视检查目标掩模 |

## 2. MASK_OUTER / MASK_INNER — 内外掩模

| 步骤 | 脚本 | 说明 |
|------|------|------|
| 2.1 | `mask_outer.py` | 外区掩模（椭圆环） |
| 2.2 | 手动编辑 | DS9 微调 |
| 2.3 | `mask_inner.py` | 内区掩模 |

## 3. CUTOUT_MASK — 合成总掩模

| 步骤 | 脚本 | 说明 |
|------|------|------|
| 3.1 | `mask_segmap_all.py` | 合并所有分割图掩模 |
| 3.2 | `mask_star.py` | 亮星掩模 |
| 3.3 | `mask_companion.py` | 伴星系掩模 |
| 3.4 | `mask_total_auto.py` | 自动合成 |
| 3.5 | `mask_manual.py` | 手动修正 |
| 3.6 | `mask_total.py` | 最终总掩模 → `TOTAL_MASK_TOTAL` |

## 4. LIMIT_DEPTH — 背景测量与极限星等

| 步骤 | 脚本 | 说明 |
|------|------|------|
| 4.1 | `profile_quarters.py` | 四象限圆形孔径测光 → `LIMIT_DEPTH_ISOTAB` |
| 4.2 | `get_bkg_rmin.py` | 自动确定背景区域起止半径 |
| 4.3 | `measure_bkg_err.py` | 测量背景误差、计算极限面亮度 |
| 4.4 | `eyeball_bkg.py` | 目视检查背景区域 |

## 5. CUTOUT_IMAGE — 图像预处理

| 步骤 | 脚本 | 说明 |
|------|------|------|
| 5.1 | `subtract_bkg.py` | 减去局部背景 → `IMG_DIR_2` |
| 5.2 | `sum_all.py` | g+r+i 合成 a 波段 |
| 5.3 | `eyeball_saturate.py` | 饱和检查 |

## 6. PROPS_ORIGINAL — 等强度线拟合与目视检查

| 步骤 | 脚本 | 说明 |
|------|------|------|
| 6.1 | `free_fitting_original.py` | 自由 Ellipse 拟合 + 矩估计 + 平滑 → `PROPS_ORIGINAL_ISOTAB` |
| 6.2 | `eyeball_fitting.py` | 1Re/2Re 目视检查图 → `PROPS_ORIGINAL_EYEBALL_FITTING` |

## 核心依赖

```
config.py          — 路径配置（所有脚本依赖）
my_tools/          — 通用工具包
  _io.py           — check_dir, run_multi
  _image.py        — 图像处理
  _plot.py         — 画图
  _cosmo.py        — 宇宙学
  _bkg.py          — 背景测量
  _fitting.py      — 拟合 + 平滑
```
