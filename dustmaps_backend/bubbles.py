from datetime import datetime, timezone
from typing import Annotated, Literal
from urllib.parse import urlsplit

from fastapi import APIRouter, Depends, HTTPException, Query
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator
from sqlalchemy import func, insert, select
from sqlalchemy.orm import Session

from .database import bubbles, operation_logs, session_dependency

router = APIRouter(prefix="/bubbles", tags=["Saved bubbles"])
DB = Annotated[Session, Depends(session_dependency)]


class BubbleInput(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)
    group: str = Field(default="default", min_length=1, max_length=50)
    lon_center: float = Field(ge=0, le=360)
    lat_center: float = Field(ge=-90, le=90)
    bubble_diameter: float = Field(gt=0, le=180)
    d_min: float = Field(ge=0)
    d_max: float = Field(gt=0)
    colormap: str | None = Field(default=None, max_length=50)
    norm: Literal["linear", "log", "hist"] | None = None
    smoothing_sigma: float | None = Field(default=None, ge=0)
    vmin: float | None = None
    vmax: float | None = None
    mark_region: bool = False
    mark_color: str | None = Field(default=None, max_length=50)
    result_image_url: str | None = Field(default=None, max_length=512)

    @field_validator("result_image_url")
    @classmethod
    def valid_url(cls, value):
        if value is not None:
            parsed = urlsplit(value)
            if not (value.startswith("/") and not value.startswith("//")) and not (parsed.scheme in {"http", "https"} and parsed.netloc):
                raise ValueError("Expected an HTTP(S) URL or an absolute URL path")
        return value

    @model_validator(mode="after")
    def ranges(self):
        if self.d_max <= self.d_min:
            raise ValueError("d_max must exceed d_min")
        if self.vmin is not None and self.vmax is not None and self.vmax <= self.vmin:
            raise ValueError("vmax must exceed vmin")
        return self


@router.post("", status_code=201)
def save_bubble(body: BubbleInput, db: DB):
    now = datetime.now(timezone.utc)
    values = body.model_dump() | {"created_at": now}
    try:
        result = db.execute(insert(bubbles).values(**values))
        bubble_id = result.inserted_primary_key[0]
        db.execute(insert(operation_logs).values(from_module="dustmaps",
            operation="[成功] 保存气泡结果", operation_time=now, operator_id=None))
        db.commit()
    except Exception:
        db.rollback()
        raise
    return {"id": bubble_id, **values}


@router.get("/groups")
def list_groups(db: DB, page: Annotated[int, Query(ge=1)] = 1,
                page_size: Annotated[int, Query(ge=1, le=100)] = 20):
    name = func.coalesce(func.nullif(bubbles.c.group, ""), "default")
    grouped = select(name.label("name"), func.count().label("count")).group_by(name).subquery()
    total = db.scalar(select(func.count()).select_from(grouped))
    rows = db.execute(select(grouped).order_by(grouped.c.name).offset((page - 1) * page_size).limit(page_size)).mappings()
    return {"items": [dict(row) for row in rows], "total": total, "page": page, "page_size": page_size}


@router.get("")
def list_bubbles(db: DB, group: str | None = None,
                 page: Annotated[int, Query(ge=1)] = 1,
                 page_size: Annotated[int, Query(ge=1, le=100)] = 20):
    query = select(bubbles)
    if group is not None:
        query = query.where(func.coalesce(func.nullif(bubbles.c.group, ""), "default") == group)
    total = db.scalar(select(func.count()).select_from(query.subquery()))
    rows = db.execute(query.order_by(bubbles.c.created_at.desc(), bubbles.c.id.desc()).offset((page - 1) * page_size).limit(page_size)).mappings()
    return {"items": [dict(row) for row in rows], "total": total, "page": page, "page_size": page_size}


@router.get("/{bubble_id:int}")
def get_bubble(bubble_id: int, db: DB):
    row = db.execute(select(bubbles).where(bubbles.c.id == bubble_id)).mappings().first()
    if row is None:
        raise HTTPException(404, "Bubble not found")
    return dict(row)
