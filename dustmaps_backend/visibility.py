"""Altitude curves in an explicit civil timezone; all plotting runs in a worker."""

import asyncio
from datetime import date, datetime, time, timedelta, timezone
from pathlib import Path
from typing import Literal
from uuid import uuid4
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

from fastapi import APIRouter, Request
from pydantic import BaseModel, ConfigDict, Field, model_validator

from .config import get_settings

router = APIRouter(tags=["visibility"])


class VisibilityInput(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False, json_schema_extra={
        "examples": [{"start_date": "2026-10-08", "end_date": "2026-10-08",
                      "timezone": "Asia/Shanghai", "observatory": "xinglong", "target": "Polaris"}]
    })
    start_date: date = Field(ge=date(1900, 1, 1), le=date(2100, 12, 30))
    end_date: date = Field(ge=date(1900, 1, 1), le=date(2100, 12, 30))
    timezone: str = Field(default="Asia/Shanghai", max_length=64, description="IANA civil timezone")
    observatory: str | None = Field(default=None, max_length=64)
    longitude: float | None = Field(default=None, ge=-180, le=180)
    latitude: float | None = Field(default=None, ge=-90, le=90)
    altitude_m: float = Field(default=0, ge=-500, le=10000)
    target: str | None = Field(default=None, max_length=64)
    coord_system: Literal["icrs", "galactic"] = "icrs"
    lon: float | None = Field(default=None, ge=0, lt=360, description="RA or Galactic longitude, degrees")
    lat: float | None = Field(default=None, ge=-90, le=90, description="Declination or Galactic latitude, degrees")
    add_moon: bool = True

    @model_validator(mode="after")
    def validate_inputs(self):
        if not 0 <= (self.end_date - self.start_date).days <= 6:
            raise ValueError("Date range must contain between one and seven inclusive dates")
        try:
            ZoneInfo(self.timezone)
        except (ZoneInfoNotFoundError, ValueError) as exc:
            raise ValueError("Unknown IANA timezone") from exc
        if self.observatory is not None:
            if self.observatory not in OBSERVATORIES:
                raise ValueError("Unknown observatory")
        elif self.longitude is None or self.latitude is None:
            raise ValueError("Select an observatory or supply longitude and latitude")
        if self.target is not None:
            if self.target not in CELESTIAL_OBJECTS:
                raise ValueError("Unknown target")
        elif self.lon is None or self.lat is None:
            raise ValueError("Select a target or supply lon and lat")
        return self


def utc_interval(params: VisibilityInput) -> tuple[datetime, datetime]:
    """Convert local midnights separately so DST days contain 23/25 hours."""
    tz = ZoneInfo(params.timezone)
    start = datetime.combine(params.start_date, time.min, tzinfo=tz)
    end = datetime.combine(params.end_date + timedelta(days=1), time.min, tzinfo=tz)
    return start.astimezone(timezone.utc), end.astimezone(timezone.utc)


def calculate_visibility(params: VisibilityInput, output_dir: str, font_path: str = "") -> dict:
    import astropy.units as u
    import matplotlib
    import numpy as np
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
    from astropy.time import Time
    from astropy.utils import iers

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.font_manager import FontProperties

    # Local tables keep the endpoint usable without network access. Predictions
    # outside the bundled IERS interval retain Astropy's degraded-accuracy warning.
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None
    site = OBSERVATORIES[params.observatory] if params.observatory else {
        "name": "Custom observatory", "lon": params.longitude,
        "lat": params.latitude, "alt": params.altitude_m}
    location = EarthLocation.from_geodetic(site["lon"] * u.deg, site["lat"] * u.deg, site["alt"] * u.m)
    start, end = utc_interval(params)
    hours = (end - start).total_seconds() / 3600
    elapsed = np.linspace(0, hours, 200 * ((params.end_date - params.start_date).days + 1))
    times = Time(start) + elapsed * u.hour
    frame = AltAz(obstime=times, location=location)
    sun_alt = get_sun(times).transform_to(frame).alt.deg
    coordinates = None
    if params.target == "Sun":
        target_alt = sun_alt
        name = CELESTIAL_OBJECTS["Sun"]["name"]
    else:
        if params.target:
            target = CELESTIAL_OBJECTS[params.target]
            sky = SkyCoord(target["ra"] * u.deg, target["dec"] * u.deg, frame="icrs")
            name = target["name"]
        else:
            sky = SkyCoord(params.lon * u.deg, params.lat * u.deg, frame=params.coord_system)
            name = "Custom target"
        target_alt = sky.transform_to(frame).alt.deg
        coordinates = {"ra_deg": float(sky.icrs.ra.deg), "dec_deg": float(sky.icrs.dec.deg),
                       "l_deg": float(sky.galactic.l.deg), "b_deg": float(sky.galactic.b.deg)}
    font = FontProperties(fname=font_path) if font_path and Path(font_path).is_file() else FontProperties()
    # English labels when no Chinese font is configured keep default installations legible.
    display_name = name if font_path else (params.target or name)
    display_site = site["name"] if font_path else (params.observatory or site["name"])
    fig, ax = plt.subplots(figsize=(8 + 2.5 * ((params.end_date - params.start_date).days + 1), 6), tight_layout=True)
    try:
        ax.fill_between(elapsed, 0, 90, where=sun_alt < 0, color="lightgrey", alpha=0.5)
        ax.fill_between(elapsed, 0, 90, where=sun_alt < -18, color="darkgrey", alpha=0.5)
        if params.add_moon:
            moon_alt = get_body("moon", times).transform_to(frame).alt.deg
            ax.plot(elapsed, moon_alt, color="gold", ls="--", label="Moon")
        ax.plot(elapsed, target_alt, lw=2, label=display_name)
        ax.set_ylim(0, 90)
        ax.set_xlim(0, hours)
        ax.set_xlabel(f"Elapsed hours since {params.start_date} 00:00 ({params.timezone})")
        ax.set_ylabel("Altitude (degrees)")
        ax.set_title(f"{params.start_date} to {params.end_date} @ {display_site}", fontproperties=font)
        ax.set_xticks(np.arange(0, hours + 1, 12))
        ax.set_xticks(np.arange(0, hours + 1, 3), minor=True)
        ax.grid(which="both", linestyle=":", linewidth=0.5)
        ax.legend(prop=font, loc="upper left", bbox_to_anchor=(1.05, 1))
        right = ax.twinx()
        ticks = ax.get_yticks()
        right.set_ylim(ax.get_ylim())
        right.set_yticks(ticks)
        with np.errstate(divide="ignore"):
            airmass = 1 / np.sin(np.deg2rad(ticks))
        right.set_yticklabels(["" if altitude <= 0 or altitude > 90 else
                              ">10" if mass > 10 else f"{mass:.2f}"
                              for altitude, mass in zip(ticks, airmass)])
        right.set_ylabel("Airmass", rotation=270, labelpad=20)
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        filename = f"visibility-{uuid4().hex}.png"
        fig.savefig(output / filename, dpi=120, bbox_inches="tight")
    finally:
        plt.close(fig)
    return {"filename": filename, "target": name, "observatory": site,
            "timezone": params.timezone, "start_utc": start.isoformat(), "end_utc": end.isoformat(),
            "coordinates": coordinates, "elapsed_hours": elapsed.tolist(),
            "altitude_deg": target_alt.tolist(), "sun_altitude_deg": sun_alt.tolist(),
            "moon_altitude_deg": moon_alt.tolist() if params.add_moon else None,
            "iers_mode": "bundled tables; accuracy degrades outside their validity interval"}


@router.get("/metadata/observatories")
def observatories():
    return {key: {**site, "timezone": "Asia/Shanghai"} for key, site in OBSERVATORIES.items()}


@router.get("/metadata/targets")
def targets():
    return CELESTIAL_OBJECTS


@router.post("/visibility/calculate")
async def visibility(params: VisibilityInput, request: Request):
    settings = get_settings()
    result = await asyncio.get_running_loop().run_in_executor(
        request.app.state.compute_pool, calculate_visibility, params,
        str(settings.output_dir), settings.font_path)
    result["plot_url"] = settings.public_base_url.rstrip("/") + "/files/" + result.pop("filename")
    return result
OBSERVATORIES = {
    'xinglong': {'name': '兴隆观测站 (河北)', 'lon': 117.58, 'lat': 40.39, 'alt': 950},
    'lijiang': {'name': '丽江/高美古观测站 (云南)', 'lon': 100.2, 'lat': 26.42, 'alt': 3200},
    'ali': {'name': '阿里观测站 (西藏)', 'lon': 80.02, 'lat': 32.32, 'alt': 5050},
    'lenghu': {'name': '冷湖天文观测基地 (青海)', 'lon': 94.27, 'lat': 38.62, 'alt': 4200},
    'nanshan': {'name': '南山观测站 (乌鲁木齐)', 'lon': 87.18, 'lat': 43.47, 'alt': 2080},
    'beijing': {'name': '北京', 'lon': 116.4, 'lat': 39.9, 'alt': 50},
    'shanghai': {'name': '上海 (佘山)', 'lon': 121.19, 'lat': 31.09, 'alt': 100},
    'guangzhou': {'name': '广州', 'lon': 113.27, 'lat': 23.13, 'alt': 10},
    'mohe': {'name': '漠河 (黑龙江)', 'lon': 122.53, 'lat': 52.97, 'alt': 500},
    'sanya': {'name': '三亚 (海南)', 'lon': 109.51, 'lat': 18.25, 'alt': 10},
}

CELESTIAL_OBJECTS = {
    'Sun': {'name': '太阳 (Sun)', 'ra': 0, 'dec': 0},
    'Polaris': {'name': '北极星 (Polaris)', 'ra': 37.95, 'dec': 89.26},
    'Sirius': {'name': '天狼星 (Sirius)', 'ra': 101.29, 'dec': -16.72},
    'AlphaCentauri': {'name': '南门二 (Alpha Centauri)', 'ra': 219.90, 'dec': -60.83},
    'M31': {'name': '仙女座星系 (M31)', 'ra': 10.68, 'dec': 41.27},
    'M42': {'name': '猎户座大星云 (M42)', 'ra': 83.82, 'dec': -5.39},
    'M45': {'name': '昴星团 (M45)', 'ra': 56.75, 'dec': 24.11},
    'M13': {'name': '武仙座球状星团 (M13)', 'ra': 250.42, 'dec': 36.46},
    'M51': {'name': '涡状星系 (M51)', 'ra': 202.47, 'dec': 47.20},
    'M8': {'name': '礁湖星云 (M8)', 'ra': 270.92, 'dec': -24.39},
    'M27': {'name': '哑铃星云 (M27)', 'ra': 299.90, 'dec': 22.72},
}