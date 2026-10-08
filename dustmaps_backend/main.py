import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from contextlib import asynccontextmanager
from datetime import datetime, timezone
from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from sqlalchemy import insert, select
from sqlalchemy.exc import SQLAlchemyError
from starlette.concurrency import run_in_threadpool
from starlette.exceptions import HTTPException as StarletteHTTPException

from . import bubbles, content
from .config import get_settings
from .database import make_engine, metadata, operation_logs

logger = logging.getLogger(__name__)
settings = get_settings()


@asynccontextmanager
async def lifespan(app: FastAPI):
    settings.output_dir.mkdir(parents=True, exist_ok=True)
    app.state.engine = make_engine(settings.database_url)
    # Spawn avoids inheriting database sockets or native scientific-library threads.
    app.state.compute_pool = ProcessPoolExecutor(
        max_workers=settings.compute_workers,
        mp_context=multiprocessing.get_context("spawn"),
    )
    app.state.active_computations = 0
    try:
        yield
    finally:
        app.state.compute_pool.shutdown(wait=True, cancel_futures=True)
        app.state.engine.dispose()


app = FastAPI(title="Dustmaps3D Backend", version="0.1.0", lifespan=lifespan,
              description="Scientific calculations and published CMS content. No user management.")
class RequestSizeLimit:
    """Count actual bytes, including requests without Content-Length."""
    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http" or scope["method"] != "POST":
            return await self.app(scope, receive, send)
        received = 0

        async def limited_receive():
            nonlocal received
            message = await receive()
            if message["type"] == "http.request":
                received += len(message.get("body", b""))
                if received > settings.max_upload_mb * 1024 * 1024:
                    raise StarletteHTTPException(413, "Request exceeds upload limit")
            return message

        await self.app(scope, limited_receive, send)


app.add_middleware(RequestSizeLimit)

if settings.cors_origins:
    app.add_middleware(CORSMiddleware, allow_origins=settings.cors_origins,
                       allow_methods=["GET", "POST"], allow_headers=["Content-Type"])


def record_operation(engine, operation, status):
    try:
        with engine.begin() as conn:
            conn.execute(insert(operation_logs).values(from_module="dustmaps",
                operation=f"[{'成功' if status < 400 else '失败'}] {operation}",
                operation_time=datetime.now(timezone.utc), operator_id=None))
    except SQLAlchemyError:
        logger.warning("Unable to persist computation operation log")


@app.middleware("http")
async def limit_computations(request: Request, call_next):
    if request.method == "POST":
        content_length = request.headers.get("content-length")
        if content_length:
            try:
                too_large = int(content_length) > settings.max_upload_mb * 1024 * 1024
            except ValueError:
                return JSONResponse({"detail": "Invalid Content-Length"}, status_code=400)
            if too_large:
                return JSONResponse({"detail": "Request exceeds upload limit"}, status_code=413)
    path = request.url.path
    compute = request.method == "POST" and (
        path.startswith("/api/v2/plots/") or path.startswith("/api/v2/dust/")
        or path.startswith("/api/v2/visibility/") or path.startswith("/api/v2/extinction/")
        or path in {"/api/v2/bubbles/analyze", "/api/v2/bubbles/schematic"})
    compute = compute or path == "/api/v2/dust/templates"
    if not compute:
        return await call_next(request)
    if request.app.state.active_computations >= settings.max_pending_computations:
        return JSONResponse({"detail": "Computation capacity reached; retry later"},
                            status_code=503, headers={"Retry-After": "5"})
    request.app.state.active_computations += 1
    status = 500
    try:
        response = await call_next(request)
        status = response.status_code
        return response
    finally:
        request.app.state.active_computations -= 1
        operations = {
            "/api/v2/dust/query": "单点计算", "/api/v2/dust/batch": "批量计算",
            "/api/v2/extinction/coefficient": "红化工具包计算",
            "/api/v2/plots/car": "绘制三维尘埃消光图", "/api/v2/plots/sin": "绘制尘埃图",
        }
        await run_in_threadpool(record_operation, request.app.state.engine, operations.get(path, path), status)


@app.exception_handler(SQLAlchemyError)
async def database_error(request: Request, exc: SQLAlchemyError):
    logger.error("Database operation failed: %s", type(exc).__name__)
    return JSONResponse({"detail": "Database unavailable or schema not initialized"}, status_code=503)


@app.get("/health", tags=["Operations"])
def health():
    return {"status": "ok", "service": "dustmaps3d-backend", "version": "0.1.0"}


@app.get("/health/ready", tags=["Operations"])
def readiness(request: Request):
    with request.app.state.engine.connect() as conn:
        for table in metadata.tables.values():
            conn.execute(select(table).limit(0))
    return {"status": "ready", "database": "connected", "datasets": {
        "dust": bool(settings.dust_data_path and Path(settings.dust_data_path).is_file()),
        "slices": bool(settings.dust_3d_path and Path(settings.dust_3d_path).is_file()),
        "dustmaps_fits": bool(settings.dustmaps_fits_path and Path(settings.dustmaps_fits_path).is_file()),
    }}


@app.get("/api/v2/viewer/config", tags=["Viewer assets"])
def viewer_config():
    if not settings.viewer_dir or not Path(settings.viewer_dir).is_dir():
        raise HTTPException(503, "Viewer data directory is not configured")
    base = settings.public_base_url.rstrip("/") + "/api/v2/viewer/assets/"
    return {"metadata_url": base + "three_volume_metadata.json",
            "bubbles_metadata_url": base + "three_volume_metadata_bubbles_abc.json",
            "brick_index_url": base + "brick_index.json", "asset_base_url": base}


@app.get("/api/v2/viewer/assets/{filename:path}", tags=["Viewer assets"])
def viewer_asset(filename: str, request: Request):
    if not settings.viewer_dir:
        raise HTTPException(503, "Viewer data directory is not configured")
    root = Path(settings.viewer_dir).resolve()
    path = (root / filename).resolve()
    if not path.is_relative_to(root) or path.suffix not in {".json", ".bin"}:
        raise HTTPException(404, "Asset not found")
    media_type = "application/json" if path.suffix == ".json" else "application/octet-stream"
    headers = {"Vary": "Accept-Encoding", "Cache-Control": "public, max-age=3600" if path.suffix == ".bin" else "no-cache"}
    if "gzip" in request.headers.get("accept-encoding", "").lower() and path.with_suffix(path.suffix + ".gz").is_file():
        path = path.with_suffix(path.suffix + ".gz")
        headers["Content-Encoding"] = "gzip"
    if not path.is_file():
        raise HTTPException(404, "Asset not found")
    return FileResponse(path, media_type=media_type, headers=headers)


app.include_router(content.router, prefix="/api/v2")
# Concrete analysis routes must precede the bubble-detail route.
from . import scientific, toolkit, visibility
app.include_router(scientific.router, prefix="/api/v2")
app.include_router(toolkit.router, prefix="/api/v2")
app.include_router(visibility.router, prefix="/api/v2")
app.include_router(bubbles.router, prefix="/api/v2")
app.mount("/files", StaticFiles(directory=settings.output_dir, check_dir=False), name="results")
if settings.cms_files_dir:
    app.mount("/cms/files", StaticFiles(directory=settings.cms_files_dir), name="cms-files")
