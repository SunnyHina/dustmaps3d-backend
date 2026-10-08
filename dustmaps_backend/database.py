"""Existing table mappings; schema creation is an explicit SQLite-only command."""
import argparse
from datetime import datetime, timezone
from pathlib import Path

from fastapi import Request

from sqlalchemy import (BigInteger, Boolean, Column, DateTime, Float, Integer,
                        MetaData, String, Table, Text, create_engine, insert, select)
from sqlalchemy.engine import make_url
from sqlalchemy.orm import Session

from .config import get_settings

metadata = MetaData()
id_type = BigInteger().with_variant(Integer, "sqlite")
categories = Table("cms_category", metadata,
    Column("id", Text, primary_key=True), Column("name", Text),
    Column("description", Text), Column("created", DateTime(timezone=True)))
articles = Table("cms_article", metadata,
    Column("id", Text, primary_key=True),
    *[Column(name, Text) for name in (
        "title", "subtitle", "abstract", "category_id", "thumbnail_path",
        "top_image_path", "content", "content_text")],
    *[Column(name, BigInteger) for name in (
        "hide_title", "show_category_article_list", "published", "ontop",
        "recommended", "deleted", "hits", "editor_id", "auditor_id")],
    *[Column(name, DateTime(timezone=True)) for name in (
        "created", "publish_time", "last_modified")])
carousels = Table("cms_carousel", metadata,
    Column("id", Text, primary_key=True),
    *[Column(name, Text) for name in ("category_id", "title", "content", "link", "image_url")],
    Column("order_", BigInteger))
bubbles = Table("saved_bubbles", metadata,
    Column("id", id_type, primary_key=True, autoincrement=True),
    Column("group", String(50)),
    *[Column(name, Float) for name in (
        "lon_center", "lat_center", "bubble_diameter", "d_min", "d_max",
        "smoothing_sigma", "vmin", "vmax")],
    Column("colormap", String(50)), Column("norm", String(20)),
    Column("mark_region", Boolean), Column("mark_color", String(50)),
    Column("result_image_url", String(512)), Column("created_at", DateTime(timezone=True)))
operation_logs = Table("operation_log", metadata,
    Column("id", id_type, primary_key=True, autoincrement=True),
    Column("from_module", Text), Column("operation", Text),
    Column("operation_time", DateTime(timezone=True)), Column("operator_id", BigInteger))


def make_engine(url: str):
    options = {"check_same_thread": False} if make_url(url).get_backend_name() == "sqlite" else {"connect_timeout": 5}
    return create_engine(url, pool_pre_ping=True, connect_args=options)


def session_dependency(request: Request):
    with Session(request.app.state.engine) as session:
        yield session


def init_sqlite(url: str, seed: bool = False):
    parsed = make_url(url)
    if parsed.get_backend_name() != "sqlite":
        raise ValueError("Initialization is SQLite-only; existing PostgreSQL schemas are never modified.")
    if parsed.database and parsed.database != ":memory:":
        Path(parsed.database).expanduser().resolve().parent.mkdir(parents=True, exist_ok=True)
    engine = make_engine(url)
    try:
        metadata.create_all(engine)
        if seed:
            now = datetime.now(timezone.utc)
            with engine.begin() as conn:
                if conn.execute(select(categories.c.id).where(categories.c.id == "demo")).first() is None:
                    conn.execute(insert(categories).values(id="demo", name="示例内容", description="SQLite 开发示例", created=now))
                if conn.execute(select(articles.c.id).where(articles.c.id == "demo-bcp")).first() is None:
                    conn.execute(insert(articles).values(id="demo-bcp", title="BCP 开发示例", category_id="demo",
                        content="<p>这是 SQLite 示例文章，不包含科学数据。</p>", content_text="SQLite 示例文章",
                        published=1, deleted=0, recommended=1, ontop=0, hits=0,
                        hide_title=0, show_category_article_list=0, created=now,
                        publish_time=now, last_modified=now))
    finally:
        engine.dispose()


def cli():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["init"])
    parser.add_argument("--seed", action="store_true", help="Add idempotent local demo content")
    args = parser.parse_args()
    try:
        init_sqlite(get_settings().database_url, args.seed)
    except ValueError as exc:
        parser.error(str(exc))
    print("SQLite initialized" + (" with demo content" if args.seed else ""))


if __name__ == "__main__":
    cli()
