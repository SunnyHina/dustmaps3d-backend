from typing import Annotated

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy import func, select
from sqlalchemy.orm import Session

from .config import get_settings
from .database import articles, carousels, categories, session_dependency

router = APIRouter(prefix="/content", tags=["Content"])
DB = Annotated[Session, Depends(session_dependency)]
Page = Annotated[int, Query(ge=1)]
PageSize = Annotated[int, Query(ge=1, le=100)]
public_article = (articles.c.published == 1) & (articles.c.deleted == 0)
article_fields = [c for c in articles.c if c.name not in ("editor_id", "auditor_id", "content_text")]
summary_fields = [c for c in article_fields if c.name != "content"]


@router.get("/categories")
def list_categories(db: DB):
    return [dict(row) for row in db.execute(select(categories).order_by(categories.c.id)).mappings()]


@router.get("/articles")
def list_articles(db: DB, category_id: str | None = None, recommended: bool = False,
                  page: Page = 1, page_size: PageSize = 20):
    condition = public_article
    if category_id is not None:
        condition &= articles.c.category_id == category_id
    if recommended:
        condition &= articles.c.recommended == 1
    total = db.scalar(select(func.count()).select_from(articles).where(condition))
    query = select(*summary_fields).where(condition).order_by(
        articles.c.ontop.desc().nulls_last(), articles.c.publish_time.desc().nulls_last(), articles.c.id.desc()
    ).offset((page - 1) * page_size).limit(page_size)
    return {"items": [dict(row) for row in db.execute(query).mappings()],
            "total": total, "page": page, "page_size": page_size}


@router.get("/bcp")
def bcp_article(db: DB):
    for article_id in get_settings().bcp_article_ids:
        row = db.execute(select(*article_fields).where(public_article, articles.c.id == article_id)).mappings().first()
        if row:
            return dict(row)
    raise HTTPException(404, "Published BCP article not found")


@router.get("/articles/{article_id}")
def article_detail(article_id: str, db: DB):
    row = db.execute(select(*article_fields).where(public_article, articles.c.id == article_id)).mappings().first()
    if row is None:
        raise HTTPException(404, "Published article not found")
    return dict(row)


@router.get("/carousels")
def list_carousels(db: DB, category_id: str | None = None):
    query = select(carousels)
    if category_id is not None:
        query = query.where(carousels.c.category_id == category_id)
    return [dict(row) for row in db.execute(query.order_by(carousels.c.order_, carousels.c.id)).mappings()]
