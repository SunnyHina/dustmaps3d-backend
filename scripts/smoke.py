"""HTTP smoke checks; never insert CMS content or saved bubbles into a shared DB."""
import argparse
import json
import time
import urllib.request


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-url", default="http://127.0.0.1:58123")
    parser.add_argument("--science", action="store_true", help="Also test configured scientific datasets and STILTS")
    args = parser.parse_args()
    base = args.base_url.rstrip("/")

    def request(path, payload=None):
        start = time.monotonic()
        data = None if payload is None else json.dumps(payload).encode()
        req = urllib.request.Request(base + path, data=data,
                                     headers={"Content-Type": "application/json"} if data else {})
        with urllib.request.urlopen(req, timeout=600) as response:
            content = response.read()
            assert response.status == 200
            result = json.loads(content) if "application/json" in response.headers.get("Content-Type", "") else content
        print(f"PASS {path} ({time.monotonic() - start:.2f}s)", flush=True)
        return result

    request("/health")
    request("/health/ready")
    request("/openapi.json")
    listing = request("/api/v2/content/articles")
    if listing["items"]:
        request("/api/v2/content/articles/" + listing["items"][0]["id"])
    request("/api/v2/content/categories")
    request("/api/v2/content/carousels")
    request("/api/v2/bubbles/groups")
    request("/api/v2/bubbles?page_size=1")
    request("/api/v2/metadata/filters")
    request("/api/v2/metadata/observatories")
    request("/api/v2/metadata/targets")
    template = request("/api/v2/dust/templates?file_format=fits")
    assert template.startswith(b"SIMPLE")
    coefficient = request("/api/v2/extinction/coefficient", {"use_2023": True, "band": "Ks"})
    assert abs(coefficient["result"] - 0.306) < 1e-10
    request("/api/v2/extinction/coefficient", {})

    def image_result(path, payload):
        result = request(path, payload)
        image_url = result.get("url", result.get("plot_url"))
        if image_url.startswith("/"):
            image_url = base + image_url
        with urllib.request.urlopen(image_url, timeout=60) as response:
            assert response.read(8) == b"\x89PNG\r\n\x1a\n"
        print("PASS PNG download", flush=True)

    image_result("/api/v2/bubbles/schematic", {"diameter": 1})
    image_result("/api/v2/visibility/calculate", {
        "start_date": "2026-10-08", "end_date": "2026-10-08", "observatory": "xinglong", "target": "Polaris"})
    if args.science:
        result = request("/api/v2/dust/query", {"coord1": 120.5, "coord2": 25.3, "d": 1.5})
        assert isinstance(result["EBV"], float)
        sky = {"d_min": 1, "d_max": 1.1, "smoothing_sigma": 0}
        image_result("/api/v2/plots/car", sky | {"lon_min": 119, "lon_max": 121, "lat_min": 24, "lat_max": 26})
        image_result("/api/v2/plots/sin", sky | {"lon_center": 120, "lat_center": 25, "fov": 2})
        image_result("/api/v2/plots/ort", {"range1_min": -0.1, "range1_max": 0.1, "range2_min": -0.1,
                                         "range2_max": 0.1, "resolution_pc": 50, "fixed_range_min": 0, "fixed_range_max": 0})
        image_result("/api/v2/plots/dpr", {})
        image_result("/api/v2/bubbles/analyze", {"l": 120.5, "b": 25.3, "d": 1, "d_low": 0.5, "d_up": 1.5, "diameter": 1})
        viewer = request("/api/v2/viewer/config")
        request(viewer["metadata_url"].removeprefix(base))


if __name__ == "__main__":
    main()
