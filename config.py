# 项目根目录
ROOT = "/data1/bjtao/StellarHalo_z02/Process2/"

# 总表
TABLE_PATH = f"{ROOT}/Catalog/proc2_limit.fits"

# cutout_mask 目录
MASK_DIR = f"{ROOT}/cutout_mask/"

# cutout_segmap 目录
SEGMAP_DIR = f"{ROOT}/cutout_segmap/"

# 各 segmap 的 mask 目录
MASK_05 = f"{SEGMAP_DIR}/segmap_05/mask"
MASK_15 = f"{SEGMAP_DIR}/segmap_15/mask"
MASK_20 = f"{SEGMAP_DIR}/segmap_20/mask"
MASK_30 = f"{SEGMAP_DIR}/segmap_30/mask"
MASK_40 = f"{SEGMAP_DIR}/segmap_40/mask"

MASK_15_TARGET = f"{SEGMAP_DIR}/segmap_15/mask_target"
MASK_30_TARGET = f"{SEGMAP_DIR}/segmap_30/mask_target"


# mask_outer 输出目录
MASK_OUTER_DIR = f"{SEGMAP_DIR}/mask_outer"
MASK_OUTER_AUTO = f"{MASK_OUTER_DIR}/mask"
MASK_OUTER_EYEBALL = f"{MASK_OUTER_DIR}/eyeball_mask"
MASK_OUTER_TARGET  = f"{MASK_OUTER_DIR}/mask_target"
MASK_OUTER_EYEBALL_TARGET = f"{MASK_OUTER_DIR}/eyeball_target"

# mask_inner 输出目录
MASK_INNER_DIR = f"{SEGMAP_DIR}/mask_inner"
MASK_INNER_AUTO = f"{MASK_INNER_DIR}/mask"
MASK_INNER_MANUAL = f"{MASK_INNER_DIR}/manual"
MASK_INNER_EYEBALL = f"{MASK_INNER_DIR}/eyeball_mask"

# mask_segmap_all 输出目录
MASK_SEGMAP_ALL_DIR = f"{MASK_DIR}/mask_segmap_all"
MASK_SEGMAP_ALL_AUTO = f"{MASK_SEGMAP_ALL_DIR}/mask"
MASK_SEGMAP_ALL_EYEBALL = f"{MASK_SEGMAP_ALL_DIR}/eyeball_mask"
MASK_SEGMAP_ALL_EYEBALL_TARGET = f"{MASK_SEGMAP_ALL_DIR}/eyeball_target"

# 图像目录
IMG_DIR = f"{ROOT}/cutout_image/remove_brightstar/"
IMG_DIR_2 = f"{ROOT}/cutout_image/subtract_background/"
VAR_DIR = f"{ROOT}/cutout_image/variance/"
EYEBALL_SATURATE_DIR = f"{ROOT}/cutout_image/eyeball_saturate/"
EYEBALL_IMG_DIR_2 = f"{ROOT}/cutout_image/subtract_background/eyeball/"

# 亮星模型与掩模输出
STAR_MODEL_DIR = f"{ROOT}/cutout_substar/star_model_all"
STAR_MASK_DIR  = f"{MASK_DIR}/mask_star/mask"
STAR_MASK_EYEBALL = f"{MASK_DIR}/mask_star/eyeball_mask"

# 伴星系掩模相关
COMPANION_CAT_DIR    = f"{ROOT}/cutout_segmap/sexcat_15/a_band"        # SExtractor 星表 (*.cat)
COMPANION_SEG_DIR    = f"{ROOT}/cutout_segmap/segmap_15/a_band"        # 分割图 (*.fits)
COMPANION_TARGET_DIR = MASK_15_TARGET  # 目标区域掩模
COMPANION_MASK_DIR   = f"{ROOT}/cutout_mask/mask_companion/mask"     # 输出目录
COMPANION_EYEBALL_DIR   = f"{ROOT}/cutout_mask/mask_companion/eyeball_mask"     # 输出目录

# 不可用的像素
UNUSABLE_MASK = f"{ROOT}/cutout_mask_big/mask_unusable"

# 总mask
TOTAL_MASK = f"{MASK_DIR}/mask_total/"
TOTAL_MASK_AUTO = f"{TOTAL_MASK}/mask_auto/"
TOTAL_EYEBALL_AUTO = f"{TOTAL_MASK}/eyeball_auto/"
TOTAL_EYEBALL_TARGET_AUTO = f"{TOTAL_MASK}/eyeball_target_auto/"
TOTAL_MASK_MANUAL = f"{TOTAL_MASK}/mask_manual/"
TOTAL_EYEBALL_MANUAL = f"{TOTAL_MASK}/eyeball_manual/"
TOTAL_REG_DS9 = f"{TOTAL_MASK}/reg_ds9/"
TOTAL_MASK_TOTAL = f"{TOTAL_MASK}/mask_total/"
TOTAL_EYEBALL_TOTAL = f"{TOTAL_MASK}/eyeball_total/"

# 估计背景误差，减背景
LIMIT_DEPTH = f"{ROOT}/limit_depth/"
LIMIT_DEPTH_ISOTAB = f"{LIMIT_DEPTH}/isotab/"
LIMIT_DEPTH_EYEBALL_BKG = f"{LIMIT_DEPTH}/eyeball_bkg/"

# original图像性质
PROPS_ORIGINAL = f"{ROOT}/props_original/"
PROPS_ORIGINAL_ISOTAB = f"{PROPS_ORIGINAL}/isotab/"
PROPS_ORIGINAL_EYEBALL_FITTING = f"{PROPS_ORIGINAL}/eyeball_fitting/"