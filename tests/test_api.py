"""Portable HTTP/database checks; never connect to the existing PostgreSQL server."""
import os
import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient
from sqlalchemy import insert, select, func

from dustmaps_backend.config import get_settings
from dustmaps_backend.database import articles, init_sqlite, make_engine, operation_logs


class APITests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        root = Path(cls.temp.name)
        cls.env = {key: os.environ.get(key) for key in ["DUSTMAPS_DATABASE_URL", "DUSTMAPS_OUTPUT_DIR", "DUSTMAPS_BCP_ARTICLE_IDS", "DUSTMAPS_MAX_UPLOAD_MB", "DUSTMAPS_DUST_DATA_PATH"]}
        os.environ.update(DUSTMAPS_DATABASE_URL=f"sqlite:///{root / 'test.db'}",
                          DUSTMAPS_OUTPUT_DIR=str(root / "results"),
                          DUSTMAPS_BCP_ARTICLE_IDS='["demo-bcp"]',
                          DUSTMAPS_MAX_UPLOAD_MB="1", DUSTMAPS_DUST_DATA_PATH="")
        get_settings.cache_clear()
        cls.url = get_settings().database_url
        init_sqlite(cls.url, seed=True)
        init_sqlite(cls.url, seed=True)
        engine = make_engine(cls.url)
        with engine.begin() as conn:
            conn.execute(insert(articles), [dict(id="draft", published=0, deleted=0),
                                            dict(id="deleted", published=1, deleted=1)])
        engine.dispose()
        from dustmaps_backend.main import app
        cls.app = app
        cls.client_context = TestClient(app)
        cls.client = cls.client_context.__enter__()

    @classmethod
    def tearDownClass(cls):
        cls.client_context.__exit__(None, None, None)
        for key, value in cls.env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
        get_settings.cache_clear()
        cls.temp.cleanup()

    def test_content_visibility_bcp_and_idempotent_initialization(self):
        self.assertEqual(self.client.get("/health/ready").status_code, 200)
        listing = self.client.get("/api/v2/content/articles").json()
        self.assertEqual(listing["total"], 1)
        self.assertNotIn("content", listing["items"][0])
        self.assertEqual(self.client.get("/api/v2/content/bcp").json()["id"], "demo-bcp")
        for identifier in ("draft", "deleted", "missing"):
            self.assertEqual(self.client.get(f"/api/v2/content/articles/{identifier}").status_code, 404)
        with self.assertRaises(ValueError):
            init_sqlite("postgresql+psycopg://invalid/database", True)

    def test_saved_bubble_transaction_and_group_pagination(self):
        body = dict(lon_center=120, lat_center=25, bubble_diameter=1, d_min=0.5, d_max=1.5)
        saved = self.client.post("/api/v2/bubbles", json=body)
        self.assertEqual(saved.status_code, 201, saved.text)
        self.assertEqual(self.client.get(f"/api/v2/bubbles/{saved.json()['id']}").status_code, 200)
        groups = self.client.get("/api/v2/bubbles/groups").json()["items"]
        self.assertEqual(groups, [{"name": "default", "count": 1}])
        self.assertEqual(self.client.get("/api/v2/bubbles?group=default").json()["total"], 1)
        self.assertEqual(self.client.post("/api/v2/bubbles", json=body | {"d_max": 0.1}).status_code, 422)
        self.assertEqual(self.client.post("/api/v2/bubbles", json=body | {"result_image_url": "javascript:alert(1)"}).status_code, 422)
        with self.app.state.engine.connect() as conn:
            self.assertEqual(conn.scalar(select(func.count()).select_from(operation_logs).where(operation_logs.c.operation == "[成功] 保存气泡结果")), 1)

    def test_streamed_payload_limit_and_busy_worker(self):
        response = self.client.post("/api/v2/dust/query", content=iter([b" " * 600000, b" " * 600000]), headers={"Content-Type": "application/json"})
        self.assertEqual(response.status_code, 413, response.text)
        multipart = iter([b'--test-boundary\r\nContent-Disposition: form-data; name="file"; filename="large.csv"\r\nContent-Type: text/csv\r\n\r\n', b"x" * 1200000, b"\r\n--test-boundary--\r\n"])
        response = self.client.post("/api/v2/dust/batch", content=multipart,
                                    headers={"Content-Type": "multipart/form-data; boundary=test-boundary"})
        self.assertEqual(response.status_code, 413, response.text)
        self.app.state.active_computations = 99
        try:
            response = self.client.get("/api/v2/dust/templates")
            self.assertEqual(response.status_code, 503)
        finally:
            self.app.state.active_computations = 0

    def test_spawn_worker_and_missing_dataset(self):
        template = self.client.get("/api/v2/dust/templates?file_format=fits")
        self.assertEqual(template.status_code, 200, template.text[:200])
        self.assertTrue(template.content.startswith(b"SIMPLE"))
        response = self.client.post("/api/v2/dust/query", json={"coord1": 120, "coord2": 25, "d": 1})
        self.assertEqual(response.status_code, 503, response.text)
        response = self.client.post("/api/v2/bubbles/schematic", json={"diameter": 1})
        self.assertEqual(response.status_code, 200, response.text)
        image = self.client.get(response.json()["url"])
        self.assertTrue(image.content.startswith(b"\x89PNG"))


if __name__ == "__main__":
    unittest.main()
