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

## 运行边界

`DUSTMAPS_COMPUTE_WORKERS` 默认 1，`DUSTMAPS_MAX_PENDING_COMPUTATIONS` 默认 4；后者限制同时处理／等待的科学计算请求及模板下载请求，超限返回 503 和 `Retry-After`。科学请求等待结果后返回，没有持久化任务队列或任务轮询接口。`DUSTMAPS_COMPUTE_TIMEOUT_SECONDS` **仅限制 STILTS 子进程**；其他计算通过输入范围、样本量、并发数量限制资源，不承诺统一执行超时。

上传默认限制 20 MiB、100000 行，可通过 `DUSTMAPS_MAX_UPLOAD_MB`、`DUSTMAPS_MAX_BATCH_ROWS` 调整。ORT 与气泡分析限制采样量。长期运行时监控科学进程内存和生成文件占用；当前不自动清理结果，尤其不能删除已保存气泡引用的图片。

安装依赖时保留 `pyproject.toml` 的版本约束：`xp-extinction-toolkit==1.0` 使用 `interp2d`，要求 `scipy<1.14`；同时限定 `numpy<2`、`pandas<3`。

`.env`、数据库、数据文件、生成文件、虚拟环境及本机部署记录不提交。代码通过标准 `pyproject.toml` 打包，uv 是开发工具而非运行时依赖。
