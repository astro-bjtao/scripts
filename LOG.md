# 工作日志

## 2026-07-01 · 管线扩展与修复

### 重构修复

**reshape_isotable bug**：重构到 `_fitting.py` 时 `iso[name].data.dtype` 改为 `data.dtype`，但 photutils 返回 `memoryview` 无 `.dtype` 属性，导致全样本 94/98 失败。修复后用 `np.array()` 统一转换 + 列级 dtype 检测。

**my_tools/__init__.py 缺失导出**：`Table`, `Column`, `fits` 在原 `my_tools.py` 顶层导入，重构时遗漏。子进程 `from my_tools import *` 拿不到 `Table`，`NameError`。已补充。

**row['survey'] str/bytes 兼容**：`run_multi_flat` 中 `str(row['survey'], encoding='utf-8')` 对 str 类型报错。改为 `isinstance(val, bytes)` 判断。

### 新增功能

| 函数 | 文件 | 说明 |
|------|------|------|
| `run_multi_flat` | `my_tools/_io.py` | 星系×波段展开调度，98×3=294任务→120进程 |

### bmodel / clean 管线整合

4 个旧脚本（各含重复函数）→ 2 个统一脚本：

| 旧 | 新 | 模式 |
|------|------|------|
| `bmodel.py` + `bmodel_var.py` | `fitting/build_model.py [image\|var]` | image→`bmodel/`, var→`bmodel_var/` |
| `clean_image.py` + `clean_var.py` | `fitting/clean_pixels.py [image\|var]` | image→`clean_image/`, var→`clean_var/` |

- 旧代码 sigma² 逻辑已移除（VAR_DIR 存的是方差非标准差）
- 新增 `config.py` 路径：`PROPS_ORIGINAL_BMODEL_VAR`, `CLEAN_IMAGE`, `CLEAN_VAR`

### 基础设施

- `compare_old/` — 重构前后对比档案（原版 my_tools.py + 55项审计脚本）
- 审计结果：所有 55 项功能测试等价
- 服务器节点选择：跑前检查 node01/02/03 负载，选最低的

### 代码规范

- **先验证后删除**：重构后必须服务器实测通过再删旧文件

### 提交（共 13 commits）

```
ef55f26 修复 run_multi_flat: row['survey'] 兼容 str/bytes
f03413c 新增 run_multi_flat: 星系×波段展开调度
6e07ed5 修复 build_model/clean_pixels: VAR_DIR 已是方差
57b7ae5 整合 bmodel/clean 管线 (4旧→2新)
3848f26 优化并行: 星系×波段展开, 跑满120核
646df8a 重构 bmodel.py: 适配当前代码框架
e67e306 修复 __init__.py: 补充 Table/Column/fits 导出
5a2aa62 存档: 重构前后 my_tools 对比档案
fdbf8ad 修复 reshape_isotable: memoryview .dtype 兼容
2441230 新增 README.md 和 LOG.md
```

## 2026-06-30 · 拟合管线自动化与代码库重构

### 拟合从手动调参到全自动化

- 4 个失败星系通过手动调参全部修复 → 提取规律
- 创建 `moments_estimate` 自动估计 eps/PA（二阶矩法）
- 创建 `smooth_isotable` 自动修复等强度线跃变和重叠
- 整合进 `free_fitting_original.py`：**98/98 零失败**

### my_tools.py 拆分 + 目录重组

- 938 行单文件 → `my_tools/` 包 7 模块
- 20 脚本平铺 → `mask/` `bkg/` `fitting/` `eyeball/` 按管线分类
- `.pth` 安装，完全向后兼容

### 关键发现

- PA 坐标系：Y+ vs X+ 差 90°
- PA 解缠：0°/180° 边界虚假跃变
- PA 漂移检测：范围 > 45° 自动扩窗

### 提交（共 12 commits）

```
2f71a9c 按管线阶段重组脚本目录结构
366cff4 拆分 my_tools.py 为功能模块
d8d9e49 PA 系统性漂移检测
1830965 PA 解缠修复
9e1ad06 smooth_isotable 等强度线平滑
a5a71ad moments_estimate 矩估计
f84cc76 整合到 free_fitting_original
```
