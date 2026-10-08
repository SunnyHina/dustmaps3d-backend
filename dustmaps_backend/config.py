from functools import lru_cache
from pathlib import Path

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    model_config = SettingsConfigDict(env_prefix="DUSTMAPS_", env_file=".env", extra="ignore")

    database_url: str = Field(default="sqlite:///./var/dustmaps.db", repr=False)
    output_dir: Path = Path("var/results")
    public_base_url: str = ""
    dust_data_path: str = ""
    dust_3d_path: str = ""
    superbubble_path: str = ""
    dustmaps_fits_path: str = ""
    font_path: str = ""
    viewer_dir: str = ""
    cms_files_dir: str = ""
    bcp_article_ids: list[str] = ["20250915131825", "20250915145027"]
    cors_origins: list[str] = []
    compute_workers: int = Field(default=1, ge=1, le=8)
    max_pending_computations: int = Field(default=4, ge=1, le=32)
    compute_timeout_seconds: int = Field(default=300, ge=1, le=3600)
    max_upload_mb: int = Field(default=20, ge=1, le=256)
    max_batch_rows: int = Field(default=100000, ge=1, le=1000000)
    stilts_command: str = "stilts"


@lru_cache
def get_settings() -> Settings:
    return Settings()
