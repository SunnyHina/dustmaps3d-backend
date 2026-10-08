# Backend facts

- Standalone FastAPI package; no Flask/chinavo_base runtime dependency. See README for install and API examples.
- Standard pyproject.toml is the source of dependency constraints; uv.lock and exported requirements.txt support uv and pip/conda environments.
- PostgreSQL uses existing CMS, saved_bubbles, operation_log tables. Startup never creates tables. SQLite initialization is explicit and rejects PostgreSQL.
- CMS and BCP are read-only published content queries; HTML is returned as data, no page rendering or user management.
- Spawned process pool handles science/plotting. Keep one HTTP worker; admission limit is per HTTP process. Scientific datasets load lazily in computation processes.
- Parquet supports point/batch/CAR/SIN; separate dustmaps3d FITS supports ORT and bubble analysis; DPR uses 3D grid and optional superbubble CSV. Scientific files and secrets stay outside Git.
- Preserve scipy<1.14 because XP extinction toolkit uses interp2d; numpy<2 and pandas<3 are also compatibility constraints.
- v2 bubble records and groups paginate separately. This deliberately differs from the legacy grouped result endpoint.
- unittest suite is self-contained; scripts/smoke.py optionally exercises configured real datasets and does not create saved bubbles or CMS content.
