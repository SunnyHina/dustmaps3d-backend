"""The original XP and 2023 extinction calculations, with explicit result units."""

import asyncio
import math
from functools import lru_cache
from typing import Literal

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel, ConfigDict, Field, model_validator

router = APIRouter(tags=["extinction"])

FILTER_GROUPS = {
    "Johnson": ["Johnson.U", "Johnson.B", "Johnson.V", "Johnson.R"],
    "Gaia": ["GAIA3.Gbp", "GAIA3.G", "GAIA3.Grp"],
    "SDSS": [f"SDSS.{b}" for b in "ugriz"],
    "Pan-STARRS 1": [f"PS1.{b}" for b in "grizy"],
    "2MASS": ["2MASS.J", "2MASS.H", "2MASS.Ks"],
    "WISE": [f"WISE.W{i}" for i in range(1, 5)],
    "Wavelength": ["440", "550"],
}
FILTER_GROUPS_2023 = {
    "Gaia": ["BP", "G", "RP"],
    "GALEX": ["FUV", "NUV"],
    "SDSS": [f"{b}'" for b in "ugriz"],
    "Pan-STARRS 1": list("grizy"),
    "2MASS": ["J", "H", "Ks"],
    "WISE": [f"W{i}" for i in range(1, 5)],
}


@lru_cache(maxsize=1)
def xp_toolkit():
    from xp_extinction_toolkit import ExtinctionToolkit

    return ExtinctionToolkit()


class ExtinctionInput(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False, json_schema_extra={
        "examples": [{"mode": "ext", "band": "GAIA3.Gbp", "ebv": 0.1, "teff": 5000},
                     {"mode": "reddening", "color": "BP-RP", "use_2023": True}]
    })
    mode: Literal["ext", "reddening"] = "ext"
    ebv: float = Field(default=0.1, ge=0, le=10, description="E(B-V), magnitudes")
    teff: float = Field(default=5000, gt=0, le=100000, description="Effective temperature, kelvin")
    band: str = Field(default="GAIA3.Gbp", max_length=64)
    color: str = Field(default="GAIA3.Gbp-GAIA3.Grp", max_length=129)
    use_2023: bool = False

    @model_validator(mode="after")
    def validate_filter(self):
        available = (set(sum(FILTER_GROUPS_2023.values(), [])) if self.use_2023
                     else set(xp_toolkit().Filters.keys()))
        identifiers = [self.band] if self.mode == "ext" else self.color.split("-")
        if len(identifiers) != (1 if self.mode == "ext" else 2) or not set(identifiers) <= available:
            raise ValueError("Unknown band/color for this model; see /metadata/filters")
        return self


def calculate_extinction(params: ExtinctionInput) -> dict:
    identifier = params.band if params.mode == "ext" else params.color
    if params.use_2023:
        from extinction_coefficient import extinction_coefficient

        result = extinction_coefficient(identifier, EBV=params.ebv, Teff=params.teff)
    else:
        tool = xp_toolkit()
        e4455 = tool.Cal_E4455(ebv=params.ebv, Teff=params.teff)
        calculate = tool.star_ext if params.mode == "ext" else tool.star_reddening
        result = calculate([identifier], E4455=e4455, Teff=[params.teff])[identifier][0]
    result = float(result)
    if not math.isfinite(result):
        raise ValueError("The scientific model returned a non-finite result")
    return {"result": result, "model": "2023" if params.use_2023 else "XP",
            "result_kind": "coefficient" if params.use_2023 else params.mode,
            "unit": "dimensionless" if params.use_2023 else "mag", "identifier": identifier}


@router.get("/metadata/filters")
def filters(use_2023: bool = False):
    groups = FILTER_GROUPS_2023 if use_2023 else FILTER_GROUPS
    available = set(sum(groups.values(), [])) if use_2023 else set(xp_toolkit().Filters.keys())
    return {"model": "2023" if use_2023 else "XP",
            "groups": [{"name": name, "filters": [b for b in bands if b in available]}
                       for name, bands in groups.items()], "filters": sorted(available),
            "color_format": "band1-band2",
            "model_limits": ("2023 clamps E(B-V) above 0.5 and temperature to each band's calibration range."
                             if use_2023 else "XP interpolates its packaged spectral grid; log(g)=4.5.")}


@router.post("/extinction/coefficient")
async def extinction(params: ExtinctionInput, request: Request):
    try:
        return await asyncio.get_running_loop().run_in_executor(
            request.app.state.compute_pool, calculate_extinction, params)
    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc
