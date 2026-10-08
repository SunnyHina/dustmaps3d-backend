"""Validated scientific APIs; CPU and plotting work runs in the process pool."""
import asyncio
import base64
import io
import math
import subprocess
import uuid
from functools import lru_cache, partial
from pathlib import Path
from typing import Annotated, Literal

from fastapi import APIRouter, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import Response
from pydantic import BaseModel, ConfigDict, Field, model_validator

from .config import get_settings

router = APIRouter(tags=["Scientific tools"])
CoordinateSystem = Literal["galactic", "equatorial"]
Longitude = Annotated[float, Field(ge=0, le=360)]
Latitude = Annotated[float, Field(ge=-90, le=90)]
Positive = Annotated[float, Field(gt=0)]
Nonnegative = Annotated[float, Field(ge=0)]


class Input(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)


class DustQuery(Input):
    coord_system: CoordinateSystem = "galactic"
    coord1: Longitude
    coord2: Latitude
    d: Positive | None = Field(default=None, description="Distance in kpc; null uses the maximum reliable distance")


class PlotStyle(Input):
    vmin: float = 0.01
    vmax: float = 2
    norm: Literal["linear", "log", "hist"] = "hist"
    colormap: Literal["Spectral_r", "viridis", "plasma", "inferno", "magma", "cividis", "Greys", "rainbow"] = "Spectral_r"

    @model_validator(mode="after")
    def colors(self):
        if self.vmin >= self.vmax or (self.norm == "log" and self.vmin <= 0):
            raise ValueError("Require vmin < vmax and a positive vmin for log normalization")
        return self


class SkyPlot(PlotStyle):
    coord_system: CoordinateSystem = "galactic"
    d_min: Nonnegative
    d_max: Positive
    smoothing_sigma: Annotated[float, Field(ge=0, le=10)] = 0.1

    @model_validator(mode="after")
    def distances(self):
        if self.d_min >= self.d_max:
            raise ValueError("d_min must be smaller than d_max")
        return self


class CarPlot(SkyPlot):
    lon_min: Longitude
    lon_max: Longitude
    lat_min: Latitude
    lat_max: Latitude

    @model_validator(mode="after")
    def bounds(self):
        if self.lat_min >= self.lat_max or self.lon_min == self.lon_max:
            raise ValueError("Require lat_min < lat_max and distinct longitude bounds")
        if self.coord_system == "equatorial" and self.lon_min > self.lon_max:
            raise ValueError("Equatorial longitude bounds must be increasing")
        return self


class SinPlot(SkyPlot):
    lon_center: Longitude
    lat_center: Latitude
    fov: Annotated[float, Field(gt=0, le=360)]
    bubble_diameter: Annotated[float, Field(gt=0, le=180)] | None = None
    mark_region: bool = False
    mark_color: Literal["black", "white", "red", "blue", "green", "yellow"] = "black"
    show_color_bar: bool = True


class OrtPlot(PlotStyle):
    axis1: Literal["x", "y", "z"] = "x"
    axis2: Literal["x", "y", "z"] = "y"
    fixed_axis: Literal["x", "y", "z"] = "z"
    range1_min: float = -2
    range1_max: float = 2
    range2_min: float = -2
    range2_max: float = 2
    fixed_range_min: float = -0.05
    fixed_range_max: float = 0.05
    resolution_pc: Positive = 10
    smooth_sigma_pc: Annotated[float, Field(ge=0, le=100)] = 5
    norm: Literal["linear", "log"] = "log"
    vmin: float = 0.01
    vmax: float = 5
    aggregate: Literal["mean", "median"] = "mean"

    @model_validator(mode="after")
    def volume(self):
        if len({self.axis1, self.axis2, self.fixed_axis}) != 3:
            raise ValueError("Axes must be a unique permutation of x, y and z")
        widths = [self.range1_max - self.range1_min, self.range2_max - self.range2_min,
                  self.fixed_range_max - self.fixed_range_min]
        step = self.resolution_pc / 1000
        if step <= 0 or not math.isfinite(step) or not all(math.isfinite(width) for width in widths):
            raise ValueError("Resolution and range widths must be finite and representable")
        if min(widths[:2]) < step or widths[2] < 0:
            raise ValueError("Plane ranges must contain at least one resolution cell; fixed range must be ordered")
        counts = [width / step for width in widths]
        if not all(math.isfinite(count) for count in counts):
            raise ValueError("Volume exceeds the supported sample count")
        points = math.prod(max(1, math.ceil(count)) for count in counts)
        if points > 2_000_000:
            raise ValueError("Volume exceeds 2,000,000 samples; increase resolution_pc or reduce ranges")
        if self.smooth_sigma_pc / self.resolution_pc > 100:
            raise ValueError("Smoothing exceeds 100 pixels; increase resolution_pc or reduce smooth_sigma_pc")
        return self


class DprPlot(PlotStyle):
    angle_degrees: Annotated[float, Field(ge=-360, le=360)] = 60
    offset: float = 0.3
    half_width: Positive = 0.1
    s_range_min: float = -3
    s_range_max: float = 3
    z_range_min: float = -0.5
    z_range_max: float = 0.5
    norm: Literal["log"] = "log"
    vmin: float = 0.05
    vmax: float = 1
    show_markers: bool = True
    show_superbubbles: bool = True

    @model_validator(mode="after")
    def bounds(self):
        if self.s_range_min >= self.s_range_max or self.z_range_min >= self.z_range_max:
            raise ValueError("Display bounds must be increasing")
        return self


class BubbleSchematic(Input):
    diameter: Annotated[float, Field(gt=0, le=180)]
    inner_factor: Annotated[float, Field(gt=0, le=1)] = 0.25
    annulus_inner_factor: Annotated[float, Field(gt=0, le=1)] = 0.375
    annulus_outer_factor: Annotated[float, Field(gt=0, le=1)] = 0.625

    @model_validator(mode="after")
    def rings(self):
        if not self.inner_factor <= self.annulus_inner_factor < self.annulus_outer_factor:
            raise ValueError("Require inner_factor <= annulus_inner_factor < annulus_outer_factor")
        if self.diameter * self.annulus_outer_factor > 90:
            raise ValueError("Outer annulus radius must not exceed 90 degrees")
        return self


class BubbleAnalysis(BubbleSchematic):
    l: Longitude
    b: Latitude
    d: Positive
    d_low: Positive
    d_up: Positive

    @model_validator(mode="after")
    def distance_range(self):
        if self.d_low >= self.d_up:
            raise ValueError("d_low must be smaller than d_up")
        radii = [math.radians(self.diameter * factor) for factor in
                 (self.inner_factor, self.annulus_inner_factor, self.annulus_outer_factor)]
        sample_estimate = 100 * 6 * 1024**2 * (1 - math.cos(radii[0]) + math.cos(radii[1]) - math.cos(radii[2]))
        if sample_estimate > 2_000_000:
            raise ValueError("Bubble region exceeds 2,000,000 samples; reduce diameter or ring factors")
        return self


@lru_cache(maxsize=3)
def load_dataset(path: str, kind: str):
    """Each compute process loads a configured dataset once, on first use."""
    import pandas as pd
    if not path or not Path(path).is_file():
        raise FileNotFoundError(f"Configure the {kind} dataset before using this endpoint")
    if kind == "dustmaps_fits":
        from astropy.table import Table
        df = Table.read(path, format="fits").to_pandas()
    else:
        df = pd.read_parquet(path)
    if kind in {"dust", "dustmaps_fits"}:
        df.index = pd.RangeIndex(len(df), name="pix_1024")
    return df


def query_fits(l, b, d, fits_path):
    """Use the package's equations with an explicit FITS dataset, without downloads."""
    import healpy as hp
    import numpy as np
    from dustmaps3d.core import read_map
    df = load_dataset(fits_path, "dustmaps_fits")
    pix = hp.ang2pix(1024, np.atleast_1d(l), np.atleast_1d(b), lonlat=True, nest=False)
    rows = df.iloc[pix].copy()
    rows["distance"] = np.atleast_1d(d)
    return read_map(rows)


def save_image(encoded, settings):
    filename = f"{uuid.uuid4().hex}.png"
    settings.output_dir.mkdir(parents=True, exist_ok=True)
    (settings.output_dir / filename).write_bytes(base64.b64decode(encoded))
    return {"filename": filename, "url": f"{settings.public_base_url.rstrip('/')}/files/{filename}"}


def execute(operation: str, params: dict, settings, upload: bytes | None = None):
    from . import science
    import numpy as np
    import pandas as pd
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy.table import Table

    settings.output_dir.mkdir(parents=True, exist_ok=True)
    query_dust = partial(query_fits, fits_path=settings.dustmaps_fits_path)
    try:
        if operation == "query":
            lon, lat = params["coord1"], params["coord2"]
            if params["coord_system"] == "equatorial":
                coordinates = SkyCoord(ra=lon * u.deg, dec=lat * u.deg, frame="icrs").galactic
                lon, lat = coordinates.l.deg, coordinates.b.deg
            results = science.calculate_dust_properties([lon], [lat], [params["d"] if params["d"] is not None else np.nan],
                                                       load_dataset(settings.dust_data_path, "dust"))
            return {"EBV" if key == "E(B-V)" else key:
                    float(values[0]) if np.isfinite(values[0]) else None for key, values in results.items()}
        if operation in {"batch", "template"}:
            if operation == "template":
                coords = ("l", "b") if params["coords"] == "galactic" else ("ra", "dec")
                df = pd.DataFrame({coords[0]: [120.5, 200.0], coords[1]: [25.3, -1.8], "d": [1.5, np.nan]})
            else:
                try:
                    if params["filename"].lower().endswith(".csv"):
                        df = pd.read_csv(io.BytesIO(upload), nrows=settings.max_batch_rows + 1)
                    else:
                        df = Table.read(io.BytesIO(upload), format="fits").to_pandas()
                except (ValueError, OSError, UnicodeError, pd.errors.ParserError) as exc:
                    raise ValueError("Cannot parse the uploaded CSV/FITS table") from exc
                if not 0 < len(df) <= settings.max_batch_rows:
                    raise ValueError(f"Upload must contain 1 to {settings.max_batch_rows} rows")
                df, error = science.validate_and_prepare_dataframe(df)
                if error:
                    raise ValueError(error)
                values = science.calculate_dust_properties(df["l"], df["b"], df["d"], load_dataset(settings.dust_data_path, "dust"))
                for source, target in {"E(B-V)": "E(B-V)_mag", "dust_density": "dust_density_mmag_pc", "sigma": "sigma_mag", "max_distance": "max_distance_kpc"}.items():
                    df[target] = values[source]
            if params["output_format"] == "csv":
                return df.to_csv(index=False).encode("utf-8")
            df = df.rename(columns={"E(B-V)_mag": "E_B_V_mag", "dust_density_mmag_pc": "dust_density", "sigma_mag": "sigma", "max_distance_kpc": "max_distance"})
            buffer = io.BytesIO()
            Table.from_pandas(df).write(buffer, format="fits")
            return buffer.getvalue()
        if operation == "car":
            result = science.plot_partial_dust_map_ic_improved(params, load_dataset(settings.dust_data_path, "dust"), settings.output_dir, settings.public_base_url)
        elif operation == "sin":
            params["radius_deg"] = params.pop("fov") / 2
            result = science.plot_partial_dust_map_ic_stilts(params, load_dataset(settings.dust_data_path, "dust"), settings.output_dir, settings.public_base_url, settings.stilts_command, settings.compute_timeout_seconds)
        elif operation == "ort":
            result = science.plot_dust_slice_web(params, query_dust)
        elif operation == "dpr":
            superbubbles = settings.superbubble_path if params.pop("show_superbubbles") else ""
            result = science.plot_oblique_slice_web(params, load_dataset(settings.dust_3d_path, "dust_3d"), superbubbles)
        elif operation == "analyze":
            return save_image(science.analyze_dust(**params, query_dust=query_dust), settings)
        elif operation == "schematic":
            return save_image(science.plot_schematic(**params), settings)
        else:
            raise RuntimeError(f"Unknown computation: {operation}")
        if result.get("error"):
            raise ValueError(result["error"])
        url = result.get("url", result.get("url_right"))
        if url.startswith("data:image/png;base64,"):
            return save_image(url.split(",", 1)[1], settings)
        return {"url": url, "filename": result["filename"]}
    finally:
        science.plt.close("all")


async def compute(request: Request, operation: str, params: dict, upload: bytes | None = None):
    try:
        return await asyncio.get_running_loop().run_in_executor(
            request.app.state.compute_pool, execute, operation, params, get_settings(), upload)
    except (FileNotFoundError, ImportError) as exc:
        raise HTTPException(503, "Required scientific dataset or dependency is unavailable; check server configuration") from exc
    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc
    except subprocess.TimeoutExpired as exc:
        raise HTTPException(504, "STILTS rendering exceeded its configured time limit") from exc


@router.post("/dust/query")
async def dust_query(body: DustQuery, request: Request):
    return await compute(request, "query", body.model_dump())


@router.post("/dust/batch")
async def dust_batch(request: Request, file: UploadFile = File(...), output_format: Literal["csv", "fits"] = Form("csv")):
    if Path(file.filename or "").suffix.lower() not in {".csv", ".fit", ".fits"}:
        raise HTTPException(422, "Only CSV and FITS files are supported")
    limit = get_settings().max_upload_mb * 1024 * 1024
    upload = await file.read(limit + 1)
    await file.close()
    if len(upload) > limit:
        raise HTTPException(413, "Upload exceeds max_upload_mb")
    content = await compute(request, "batch", {"filename": file.filename, "output_format": output_format}, upload)
    return download(content, output_format, "dust_results")


def download(content: bytes, output_format: str, name: str):
    return Response(content, media_type="text/csv" if output_format == "csv" else "application/fits",
                    headers={"Content-Disposition": f'attachment; filename="{name}.{output_format}"'})


@router.get("/dust/templates")
async def dust_template(request: Request, coords: CoordinateSystem = "galactic", file_format: Literal["csv", "fits"] = "csv"):
    content = await compute(request, "template", {"coords": coords, "output_format": file_format})
    return download(content, file_format, f"template_{coords}")


@router.post("/plots/car")
async def plot_car(body: CarPlot, request: Request):
    return await compute(request, "car", body.model_dump())


@router.post("/plots/sin")
async def plot_sin(body: SinPlot, request: Request):
    return await compute(request, "sin", body.model_dump())


@router.post("/plots/ort")
async def plot_ort(body: OrtPlot, request: Request):
    return await compute(request, "ort", body.model_dump())


@router.post("/plots/dpr")
async def plot_dpr(body: DprPlot, request: Request):
    return await compute(request, "dpr", body.model_dump())


@router.post("/bubbles/analyze")
async def bubble_analyze(body: BubbleAnalysis, request: Request):
    return await compute(request, "analyze", body.model_dump())


@router.post("/bubbles/schematic")
async def bubble_schematic(body: BubbleSchematic, request: Request):
    return await compute(request, "schematic", body.model_dump())
