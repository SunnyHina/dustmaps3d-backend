# dustmaps3d-backend

独立 FastAPI 后端，提供尘埃查询、批量计算、四类绘图、气泡分析与保存、消光／红化计算、天体可见性以及已发布 CMS 内容查询。支持 PostgreSQL 与 SQLite，通过 HTTP API 提供服务，不包含前端或用户认证模块。

`GET /api/v2/content/bcp` 返回配置指定的已发布 BCP 文章。文章 HTML 作为 JSON 字段返回，供调用方展示。

## 安装

要求 Python **3.11 或 3.12**。以下命令均在项目根目录执行。科学依赖较多，建议使用 Python 3.11；不需要 Node.js 或 npm。

### uv

```bash
uv sync --locked
cp .env.example .env
uv run dustmaps-db init --seed
```

`uv.lock` 固定依赖版本，`uv sync` 创建项目 `.venv` 并安装项目。可参考 [uv 环境同步文档](https://docs.astral.sh/uv/concepts/projects/sync/)。

### Python venv / pip

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip install --no-deps -e .
cp .env.example .env
python -m dustmaps_backend.database init --seed
```

`requirements.txt` 从锁文件导出，用于安装固定版本的依赖；`pyproject.toml` 定义项目及其依赖范围。

### conda

```bash
conda create -n dustmaps3d-backend python=3.11 pip
conda activate dustmaps3d-backend
python -m pip install -r requirements.txt
python -m pip install --no-deps -e .
cp .env.example .env
python -m dustmaps_backend.database init --seed
```

conda 管理 Python 环境，pip 安装项目及同一组 Python 依赖；参见 [conda 环境文档](https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html)。上述安装方式任选一种。已有 `.env` 时保留并修改它，不要重复覆盖。

## 数据库

`.env.example` 默认使用本地 SQLite：

```dotenv
DUSTMAPS_DATABASE_URL=sqlite:///./var/dustmaps.db
DUSTMAPS_BCP_ARTICLE_IDS=["demo-bcp"]
```

`dustmaps-db init` 显式创建所需 SQLite 表；`--seed` 幂等插入示例分类与 ID 为 `demo-bcp` 的已发布文章。初始化不会生成科学数据。服务启动不会自动建表；未初始化时数据库接口返回 503。SQLite 适合本地开发和独立部署。

使用 PostgreSQL 时，在 `.env` 中设置连接信息：

```dotenv
DUSTMAPS_DATABASE_URL=postgresql+psycopg://user:password@localhost:5432/dustmaps
DUSTMAPS_BCP_ARTICLE_IDS=["your-published-article-id"]
```

替换占位凭据及文章 ID，密码中的特殊字符需要 URL 编码。初始化命令仅适用于 SQLite；PostgreSQL 部署需要预先准备下列数据表，字段定义见 `dustmaps_backend/database.py`：

| 表 | 本项目用途 |
| --- | --- |
| `cms_category`、`cms_article`、`cms_carousel` | 只读内容；文章限定 `published=1` 且 `deleted=0` |
| `saved_bubbles` | 保存气泡参数、结果图片 URL；查询记录与分组 |
| `operation_log` | 追加计算与保存操作记录，`operator_id` 为空 |

文章 ID 是字符串。BCP 按 `DUSTMAPS_BCP_ARTICLE_IDS` 顺序选择第一篇有效文章。数据库连接与科学文件位置分别配置，连接数据库不会获得科学数据或图片。

## 科学数据与配置

配置从进程环境和根目录 `.env` 读取；进程环境优先。相对文件路径以启动时工作目录为基准。科学文件不提交到仓库，不在请求过程中自动下载；需要部署者按下述格式准备数据。

| 配置变量 | 数据或用途 |
| --- | --- |
| `DUSTMAPS_DUST_DATA_PATH` | 单点、批量、CAR、SIN 使用的模型参数 Parquet |
| `DUSTMAPS_DUSTMAPS_FITS_PATH` | ORT 与气泡密度分析使用的 `dustmaps3d` 兼容全天天区模型 FITS |
| `DUSTMAPS_DUST_3D_PATH` | DPR 使用的三维网格 Parquet，至少包含 `X`、`Y`、`Z`、`dust` |
| `DUSTMAPS_SUPERBUBBLE_PATH` | 可选超级气泡参数 CSV，供 DPR 叠加截面 |
| `DUSTMAPS_FONT_PATH` | 可选中文字体文件，供可见性图使用；留空使用英文标签 |
| `DUSTMAPS_VIEWER_DIR` | 三维 JSON 元数据、BIN 体数据及分块目录 |
| `DUSTMAPS_CMS_FILES_DIR` | 可选公开 CMS 附件目录，映射 `/cms/files/` |
| `DUSTMAPS_OUTPUT_DIR` | 生成图片目录，默认 `var/results`，映射 `/files/` |
| `DUSTMAPS_PUBLIC_BASE_URL` | 图片与三维资源 URL 的服务地址前缀，留空返回绝对 URL 路径 |
| `DUSTMAPS_CORS_ORIGINS` | 允许的前端来源 JSON 数组，例如 `["http://localhost:5173"]` |

例如：

```dotenv
DUSTMAPS_DUST_DATA_PATH=/srv/dustmaps-data/data.parquet
DUSTMAPS_DUSTMAPS_FITS_PATH=/srv/dustmaps-data/dust_model.fits
DUSTMAPS_DUST_3D_PATH=/srv/dustmaps-data/grid3d.parquet
DUSTMAPS_SUPERBUBBLE_PATH=/srv/dustmaps-data/superbubbles.csv
DUSTMAPS_VIEWER_DIR=/srv/dustmaps-data/volume
DUSTMAPS_OUTPUT_DIR=var/results
DUSTMAPS_PUBLIC_BASE_URL=http://localhost:58123
```

模型 Parquet 与 FITS 的行顺序必须对应 **NSIDE=1024、RING 排序** 的 HEALPix 像素，不能用任意截取后重新编号的表替代。Parquet 需要 `b_lim`、`bubble`、`diffuse_dust_rho`、`h`、`max_distance`、`sigma`，以及四组 `distance_1…4`、`span_1…4`、`Cum_EBV_1…4`。保留完整原始文件及字段，CAR/SIN 从 HEALPix 像素编号恢复天区坐标。DPR 网格坐标使用 kpc，`dust` 为尘埃密度；超级气泡 CSV 需要形状、中心、半轴及角度字段，读取定义见 `science.py` 中 `load_dpr_superbubbles`。

SIN 绘图还需要 Java 与可执行的 [STILTS](https://www.star.bristol.ac.uk/~mbt/stilts/) 命令。`DUSTMAPS_STILTS_COMMAND` 默认 `stilts`，可以改成可执行脚本绝对路径；不能填入带参数的整条 shell 命令。

未配置数据时，内容、气泡保存、消光计算、可见性、气泡示意图及模板下载仍可使用；依赖缺失数据的科学请求返回 503。科学数据在每个计算进程首次使用时加载并缓存，增加进程数量会增加内存占用。

## 启动服务

先确认示例高位端口 `58123` 未被占用：

```bash
ss -ltn 'sport = :58123'
mkdir -p var/matplotlib
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MPLCONFIGDIR="$PWD/var/matplotlib"
uv run uvicorn dustmaps_backend.main:app --host 127.0.0.1 --port 58123 --workers 1
```

pip／conda 环境激活后，把 `uv run uvicorn` 改成 `python -m uvicorn`。需要其他设备访问时将监听地址改为 `0.0.0.0`，按部署环境配置访问范围；项目不提供用户认证。服务应保持 **1 个 HTTP worker**，科学计算进程数量另由配置控制，避免重复加载数据及分散并发限制。

- API 文档：`http://localhost:58123/docs`；契约：`/openapi.json`。
- `GET /health` 检查进程；`GET /health/ready` 检查数据库表并报告数据文件存在状态，数据文件缺失不影响数据库就绪判断。
- 以下 SQLite 内容请求不需要科学数据：

```bash
curl --fail http://localhost:58123/health/ready
curl --fail http://localhost:58123/api/v2/content/bcp
```

## API 与输入示例

以下路径均以 `/api/v2` 开头。POST 默认 JSON，只有批量上传使用 `multipart/form-data`。完整可选字段与约束以 `/docs` 为准。

| 方法与路径 | 说明 |
| --- | --- |
| `GET /content/categories` | 分类 |
| `GET /content/articles?category_id=demo&page=1&page_size=20` | 文章摘要；支持 `recommended=true` |
| `GET /content/articles/{id}`、`GET /content/bcp` | 文章正文、BCP 文章 |
| `GET /content/carousels?category_id=demo` | 轮播 |
| `GET /metadata/filters?use_2023=false` | 对应模型的滤光片、颜色命名和边界说明 |
| `GET /metadata/observatories`、`GET /metadata/targets` | 可见性预设台站与天体 |
| `POST /dust/query`、`POST /dust/batch` | 单点查询、CSV/FITS 批量计算 |
| `GET /dust/templates?coords=galactic&file_format=csv` | 下载模板，亦支持 `equatorial` 与 `fits` |
| `POST /plots/car`、`/plots/sin`、`/plots/ort`、`/plots/dpr` | 天区图和切片 |
| `POST /bubbles/analyze`、`POST /bubbles/schematic` | 密度／峰值分析图、内圆外环示意图 |
| `POST /bubbles`、`GET /bubbles`、`GET /bubbles/{id}` | 保存、分页读取与详情 |
| `GET /bubbles/groups` | 分页查询分组及记录数量 |
| `POST /extinction/coefficient`、`POST /visibility/calculate` | 消光／红化、天体可见性 |
| `GET /viewer/config`、`GET /viewer/assets/{filename}` | 三维资源地址与 JSON/BIN 文件，支持已有 `.gz` 资源协商 |

科学 POST 示例：

| 路径 | JSON 请求体 |
| --- | --- |
| `/dust/query` | `{"coord_system":"galactic","coord1":120.5,"coord2":25.3,"d":1.5}` |
| `/plots/car` | `{"coord_system":"galactic","lon_min":120,"lon_max":121,"lat_min":25,"lat_max":26,"d_min":0.1,"d_max":1.5}` |
| `/plots/sin` | `{"lon_center":120.5,"lat_center":25.3,"fov":2,"d_min":0.1,"d_max":1.5}` |
| `/plots/ort` | `{"axis1":"x","axis2":"y","fixed_axis":"z","range1_min":-0.1,"range1_max":0.1,"range2_min":-0.1,"range2_max":0.1,"fixed_range_min":-0.01,"fixed_range_max":0.01,"resolution_pc":20}` |
| `/plots/dpr` | `{"angle_degrees":60,"offset":0.3,"half_width":0.1,"s_range_min":-3,"s_range_max":3,"z_range_min":-0.5,"z_range_max":0.5,"show_superbubbles":false}` |
| `/bubbles/analyze` | `{"l":120.5,"b":25.3,"d":1,"d_low":0.8,"d_up":1.2,"diameter":1}` |
| `/bubbles/schematic` | `{"diameter":1,"inner_factor":0.25,"annulus_inner_factor":0.375,"annulus_outer_factor":0.625}` |
| `/extinction/coefficient`，XP 消光 | `{"mode":"ext","band":"GAIA3.Gbp","ebv":0.1,"teff":5000}` |
| `/extinction/coefficient`，2023 红化系数 | `{"use_2023":true,"mode":"reddening","color":"BP-RP","ebv":0.1,"teff":5000}` |
| `/visibility/calculate` | `{"start_date":"2026-10-08","end_date":"2026-10-08","timezone":"Asia/Shanghai","observatory":"xinglong","target":"Polaris"}` |

例如发起 JSON 请求：

```bash
curl --fail http://localhost:58123/api/v2/bubbles/schematic \
  -H 'Content-Type: application/json' -d '{"diameter":1}'
```

批量文件要求 `l,b,d` 或 `ra,dec,d` 三列，列名不区分大小写；必须保留 `d` 列，空距离使用该方向最大可靠距离：

```bash
curl --fail 'http://localhost:58123/api/v2/dust/templates?coords=galactic&file_format=csv' -o input.csv
curl --fail http://localhost:58123/api/v2/dust/batch \
  -F file=@input.csv -F output_format=csv -o result.csv
```

保存绘图结果时，先读取绘图响应的 `url`，再提交实际 URL，例如：

```json
{"group":"example","lon_center":120.5,"lat_center":25.3,"bubble_diameter":1,"d_min":0.8,"d_max":1.2,"result_image_url":"/files/replace-with-generated-filename.png"}
```

`GET /bubbles?group=example&page=1&page_size=20` 按**记录**分页；`GET /bubbles/groups` 单独按**分组**分页。

单位与响应约定：

- 天区角度为度，尘埃距离、切片范围和偏移为 kpc；`resolution_pc`、`smooth_sigma_pc` 明确使用 pc。气泡 `diameter` 为角直径，环半径等于直径乘对应 factor。
- 单点响应 `EBV`、`sigma` 为 mag，`dust_density` 为 mmag/pc（数值等于 mag/kpc），`max_distance` 为 kpc；无有效值时为 `null`。
- XP 返回消光／红化量，单位 mag；2023 返回无量纲系数，不能混用。2023 超过校准范围时采用上游库边界值；滤光片名称也不同，例如 `GAIA3.Gbp` 与 `BP`。
- 绘图与分析返回 `url`、`filename`；可见性返回 `plot_url`、时间序列及高度角。生成 PNG 经 `/files/` 访问；数据库只保存 URL，不保存图像字节。
- 可见性日期包含首尾，最多七天；自定义台站用 `longitude`、`latitude`、`altitude_m`，自定义目标用 `coord_system: "icrs"` 或 `"galactic"` 和 `lon`、`lat`。时区使用 IANA 名称，正确处理夏令时；X 轴是从本地零点开始的实际经过小时。Astropy 使用内置 IERS 表，超出有效期会降低精度，响应包含说明。
- CMS 附件与数据库中的图片 URL 不自动重写。部署时应确保这些 URL 可访问，或由调用方映射到配置的附件服务。

## 绘图接口

以下 POST 接口接收 JSON，完成绘图后返回 PNG 地址。路径均以 `/api/v2` 开头；可直接使用上一节的请求体示例。

| 接口 | 图片内容 | 必要数据／依赖 |
| --- | --- | --- |
| `POST /plots/car` | 指定经纬度范围、距离区间的尘埃天区图 | 模型 Parquet |
| `POST /plots/sin` | 以中心点和视场选择天区，支持气泡区域标记 | 模型 Parquet、Java、STILTS |
| `POST /plots/ort` | 银河笛卡尔坐标系中的正交切片 | 模型 FITS |
| `POST /plots/dpr` | 给定角度、偏移和厚度的斜切片，可叠加超级气泡截面 | 三维网格 Parquet；叠加时配置气泡 CSV |
| `POST /bubbles/analyze` | 气泡内圆与外环的密度曲线、差值及峰值分析图 | 模型 FITS |
| `POST /bubbles/schematic` | 内圆和外环的几何示意图 | 无科学数据依赖 |
| `POST /visibility/calculate` | 目标、太阳及可选月球的高度角随时间变化图 | 无尘埃数据依赖 |

### 天区图：CAR 与 SIN

两种天区图都必须提供 `d_min`、`d_max`，单位 kpc，满足 `0 <= d_min < d_max`。`coord_system` 默认为 `galactic`（银道坐标），也可选 `equatorial`（赤道坐标）。经度范围为 `[0, 360]`，纬度范围为 `[-90, 90]`，角度单位均为度。

| 接口 | 参数 | 含义与约束 |
| --- | --- | --- |
| CAR | `lon_min`、`lon_max`、`lat_min`、`lat_max`，必填 | 经纬度边界；纬度必须递增，经度不能相同。银道经度支持跨零点，如 `350 → 10`；赤道经度必须递增 |
| SIN | `lon_center`、`lat_center`、`fov`，必填 | 中心坐标与视场直径；`0 < fov <= 360`，选区半径为 `fov / 2` |
| SIN | `mark_region=false`、`bubble_diameter=null` | 启用区域标记；指定气泡角直径时以其一半为标记半径，否则使用 `fov / 5`。角直径范围 `(0, 180]` |
| SIN | `mark_color="black"`、`show_color_bar=true` | 标记颜色可选 `black`、`white`、`red`、`blue`、`green`、`yellow`；控制色条显示 |
| 两者 | `smoothing_sigma=0.1` | 角平滑尺度，范围 `[0, 10]` 度，`0` 表示不平滑 |

SIN **只接收中心点与视场，不接收经纬度边界字段**；调用方需要自行确定中心和视场。CAR 与 SIN 的投影不同。SIN 在 `fov >= 180` 且未指定 `bubble_diameter` 时使用全天 Aitoff 投影；局部图会在选区半径外增加显示余量。

### 正交切片：ORT

`POST /plots/ort` 的所有参数均有默认值，位置与范围使用 kpc：

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `axis1`、`axis2`、`fixed_axis` | `"x"`、`"y"`、`"z"` | 横轴、纵轴、厚度方向；必须是 `x/y/z` 的不重复排列 |
| `range1_min/max`、`range2_min/max` | 各 `-2 / 2` | 绘图平面两个方向的范围，必须递增 |
| `fixed_range_min/max` | `-0.05 / 0.05` | 沿第三轴聚合的范围；两端相等表示单平面 |
| `resolution_pc` | `10` | 采样间隔，单位 **pc**，必须大于零；越小采样越密 |
| `smooth_sigma_pc` | `5` | 平滑尺度，单位 **pc**，范围 `[0, 100]` |
| `aggregate` | `"mean"` | 沿厚度方向使用均值 `mean` 或中位数 `median` |

每个平面方向至少包含一个采样间隔，估算采样总数不得超过 2,000,000；超限时增大 `resolution_pc` 或缩小范围。`smooth_sigma_pc / resolution_pc` 不得超过 100。

### 斜切片：DPR

`POST /plots/dpr` 的所有参数均有默认值：

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `angle_degrees` | `60` | XY 平面内切片方向角，范围 `[-360, 360]` 度；沿平面方向为 `(cos θ, sin θ)` |
| `offset` | `0.3` | 沿法向 `(-sin θ, cos θ)` 的有符号偏移，单位 kpc |
| `half_width` | `0.1` | 切片半厚度，单位 kpc，必须大于零；完整厚度为其两倍 |
| `s_range_min/max` | `-3 / 3` | 沿切片方向的显示范围，单位 kpc |
| `z_range_min/max` | `-0.5 / 0.5` | Z 方向显示范围，单位 kpc |
| `show_markers` | `true` | 显示图中 `(s=0, z=0)` 参考标记 |
| `show_superbubbles` | `true` | 叠加配置 CSV 中的气泡截面；设为 `false` 可关闭 |

两个显示范围均须递增。DPR 使用固定 15 pc 网格，不接收 ORT 的 `resolution_pc` 或 `smooth_sigma_pc` 参数；显示范围不改变网格分辨率。未配置气泡 CSV 时不绘制气泡叠加层。

### 尘埃图颜色设置

以下参数适用于 CAR、SIN、ORT、DPR，不适用于气泡分析／示意图和可见性图。

| 参数 | CAR / SIN | ORT | DPR |
| --- | --- | --- | --- |
| `norm` | 默认 `hist`；可选 `linear`、`log`、`hist` | 默认 `log`；可选 `linear`、`log` | 仅 `log` |
| `vmin` / `vmax` | `0.01 / 2` | `0.01 / 5` | `0.05 / 1` |
| `colormap` | 默认 `Spectral_r` | 默认 `Spectral_r` | 默认 `Spectral_r` |

`vmin`、`vmax` 控制密度色标范围，须满足 `vmin < vmax`，对数色标还要求 `vmin > 0`。`colormap` 可选 `Spectral_r`、`viridis`、`plasma`、`inferno`、`magma`、`cividis`、`Greys`、`rainbow`。

### 气泡分析与示意图

两个接口共用 `diameter`（必填，角直径，单位度，范围 `(0, 180]`）及三个半径系数：`inner_factor=0.25`、`annulus_inner_factor=0.375`、`annulus_outer_factor=0.625`。半径等于 `diameter × factor`，并满足 `0 < inner_factor <= annulus_inner_factor < annulus_outer_factor <= 1`，外环外半径不超过 90 度。

`POST /bubbles/schematic` 仅需要上述几何参数。`POST /bubbles/analyze` 还必须提供：

| 参数 | 含义 |
| --- | --- |
| `l`、`b` | 银经、银纬，单位度；仅接受银道坐标 |
| `d` | 正数参考距离，单位 kpc；当前为必填字段，但不参与分析计算 |
| `d_low`、`d_up` | 峰值分析的距离窗口，单位 kpc，满足 `0 < d_low < d_up` |

分析图展示内圆／外环平均密度及差值，并标注距离窗口中的峰值结果；响应只提供图片地址，不包含独立的数值拟合结果。分析采样量限制为 2,000,000，超限时缩小角直径或区域系数。绘图不会自动创建气泡记录；保存时需另行调用 `POST /bubbles`。

### 天体可见性图

`POST /visibility/calculate` 必填 `start_date`、`end_date`，格式 `YYYY-MM-DD`，包含首尾日期，最多七天；支持日期范围 `1900-01-01` 至 `2100-12-30`。`timezone` 默认 `Asia/Shanghai`，使用 IANA 时区名。

- 台站：通过 `GET /metadata/observatories` 获取 `observatory` 预设，或提供 `longitude`、`latitude`（度）和 `altitude_m`（米，默认 `0`）。
- 目标：通过 `GET /metadata/targets` 获取 `target` 预设，或提供 `coord_system`（`icrs` 或 `galactic`）、`lon`、`lat`（度）。
- `add_moon=true` 默认加入月球曲线；设为 `false` 时不计算月球高度角。

响应含 `plot_url`、`elapsed_hours`、`altitude_deg`、`sun_altitude_deg`、`moon_altitude_deg` 等字段，可直接展示图片或使用序列自行绘制；关闭月球时 `moon_altitude_deg` 为 `null`。计算使用内置 IERS 表，超出其有效期时精度会降低。

### 请求、返回与图片访问

以 SIN 绘图为例：

```bash
curl --fail http://localhost:58123/api/v2/plots/sin \
  -H 'Content-Type: application/json' \
  -d '{"lon_center":120.5,"lat_center":25.3,"fov":2,"d_min":0.1,"d_max":1.5,"mark_region":true,"bubble_diameter":1}'
```

四类尘埃图及两类气泡图返回相同结构（文件名仅为示例）：

```json
{"url":"/files/generated-image.png","filename":"generated-image.png"}
```

随后向响应中的 `url` 发起 GET 请求即可读取或下载 PNG。可见性图改用 `plot_url` 字段。POST 本身返回 JSON，不返回图片字节或 Base64。未设置 `DUSTMAPS_PUBLIC_BASE_URL` 时，`/files/...` 应相对于**后端服务地址**解析；前后端域名不同时，配置该变量或由前端拼接后端地址。结果文件目录由 `DUSTMAPS_OUTPUT_DIR` 决定。

所有绘图请求均拒绝未定义字段。参数错误或超过采样限制返回 `422`；科学数据／依赖不可用或计算并发达到上限返回 `503`；STILTS 超时返回 `504`。详细错误见响应的 `detail` 字段；并发超限时同时返回 `Retry-After`。

### 三维交互展示

`GET /viewer/config` 提供元数据、气泡元数据、分块索引及资源地址；`GET /viewer/assets/{filename}` 提供 JSON／BIN 文件。部署时配置 `DUSTMAPS_VIEWER_DIR`。这些接口不生成三维视图 PNG；体渲染、气泡三维展示、旋转与缩放由前端实现。

## 运行边界

`DUSTMAPS_COMPUTE_WORKERS` 默认 1，`DUSTMAPS_MAX_PENDING_COMPUTATIONS` 默认 4；后者限制同时处理／等待的科学计算请求及模板下载请求，超限返回 503 和 `Retry-After`。科学请求等待结果后返回，没有持久化任务队列或任务轮询接口。`DUSTMAPS_COMPUTE_TIMEOUT_SECONDS` **仅限制 STILTS 子进程**；其他计算通过输入范围、样本量、并发数量限制资源，不承诺统一执行超时。

上传默认限制 20 MiB、100000 行，可通过 `DUSTMAPS_MAX_UPLOAD_MB`、`DUSTMAPS_MAX_BATCH_ROWS` 调整。ORT 与气泡分析限制采样量。长期运行时监控科学进程内存和生成文件占用；当前不自动清理结果，尤其不能删除已保存气泡引用的图片。

安装依赖时保留 `pyproject.toml` 的版本约束：`xp-extinction-toolkit==1.0` 使用 `interp2d`，要求 `scipy<1.14`；同时限定 `numpy<2`、`pandas<3`。

`.env`、数据库、数据文件、生成文件、虚拟环境及本机部署记录不提交。代码通过标准 `pyproject.toml` 打包，uv 是开发工具而非运行时依赖。
