import io
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import healpy as hp
import numpy as np
import pandas as pd
from astropy.table import Table
from pydantic import ValidationError

from dustmaps_backend import science
from dustmaps_backend.scientific import BubbleSchematic, CarPlot, DustQuery, OrtPlot, SinPlot, execute


class ScientificTests(unittest.TestCase):
    def test_sin_accepts_negative_longitude_and_rejects_nonfinite_centers(self):
        base = dict(lat_center=0, fov=10, d_min=0.3, d_max=0.5)
        for frame in ("galactic", "equatorial"):
            for longitude in (-180, -5, 0, 355, 360):
                with self.subTest(frame=frame, longitude=longitude):
                    plot = SinPlot(coord_system=frame, lon_center=longitude, **base)
                    self.assertEqual(plot.lon_center, longitude)
        for longitude in (float("nan"), float("inf"), float("-inf")):
            with self.subTest(longitude=longitude), self.assertRaises(ValidationError):
                SinPlot(lon_center=longitude, **base)

    def test_car_longitude_bounds_preserve_negative_and_full_sky_ranges(self):
        base = dict(lat_min=-5, lat_max=5, d_min=0.3, d_max=0.5)
        for frame in ("galactic", "equatorial"):
            for low, high in [(-5, 5), (-15, -5), (-180, 180), (0, 360)]:
                with self.subTest(frame=frame, low=low, high=high):
                    plot = CarPlot(coord_system=frame, lon_min=low, lon_max=high, **base)
                    self.assertEqual((plot.lon_min, plot.lon_max), (low, high))
        self.assertEqual(CarPlot(lon_min=355, lon_max=5, **base).lon_min, 355)
        for bounds in [dict(lon_min=float("nan"), lon_max=5),
                       dict(lon_min=-5, lon_max=float("inf")),
                       dict(lon_min=-5, lon_max=-5),
                       dict(coord_system="equatorial", lon_min=355, lon_max=5)]:
            with self.subTest(bounds=bounds), self.assertRaises(ValidationError):
                CarPlot(**bounds, **base)

    def test_model_derivative_and_blank_distance_batch_round_trip(self):
        # A synthetic, nonzero four-cloud sightline checks the numerical equations
        # independently of the production dataset and upload serialization.
        row = dict(b_lim=20.0, bubble=0.02, diffuse_dust_rho=0.5, h=0.1,
                   max_distance=2.0, sigma=0.04)
        for n in range(1, 5):
            row.update({f"distance_{n}": n * 0.3, f"span_{n}": 0.1, f"Cum_EBV_{n}": 0.2})
        params = [row[k] for k in ["b_lim", "bubble", "diffuse_dust_rho", "h"]]
        params += [row[f"{key}_{n}"] for n in range(1, 5) for key in ["distance", "span", "Cum_EBV"]]
        step, distance = 1e-6, 1.1
        derivative = (science.component4(distance + step, *params) - science.component4(distance - step, *params)) / (2 * step)
        self.assertAlmostEqual(float(derivative), float(science.derivative_of_component4(distance, *params)), places=6)
        pixel = hp.ang2pix(1024, 120.5, 25.3, lonlat=True)
        dataset = pd.DataFrame([row], index=[pixel])
        with tempfile.TemporaryDirectory() as directory:
            settings = SimpleNamespace(output_dir=Path(directory), dustmaps_fits_path="", dust_data_path="synthetic", max_batch_rows=3)
            with patch("dustmaps_backend.scientific.load_dataset", return_value=dataset):
                encoded = execute("batch", {"filename": "input.csv", "output_format": "fits"}, settings, b"l,b,d\n120.5,25.3,\n120.5,25.3,1.1\n")
                table = Table.read(io.BytesIO(encoded), format="fits").to_pandas()
                self.assertAlmostEqual(table["E_B_V_mag"][0], float(science.component4(2.0, *params)))
                self.assertAlmostEqual(table["dust_density"][1], float(derivative), places=6)
                self.assertEqual(table["max_distance"].tolist(), [2.0, 2.0])
                single = execute("query", {"coord_system": "galactic", "coord1": 120.5, "coord2": 25.3, "d": None}, settings)
                self.assertEqual(single["EBV"], table["E_B_V_mag"][0])

    def test_resource_and_input_validation(self):
        for constructor, data in [(DustQuery, {"coord1": 361, "coord2": 0}),
                                  (DustQuery, {"coord1": 10, "coord2": 0, "d": float("inf")}),
                                  (OrtPlot, {"resolution_pc": 0.001}),
                                  (OrtPlot, {"resolution_pc": 5e-324}),
                                  (OrtPlot, {"range1_min": -1e308, "range1_max": 1e308}),
                                  (OrtPlot, {"resolution_pc": 1e-300}),
                                  (OrtPlot, {"axis2": "x"}),
                                  (BubbleSchematic, {"diameter": 1, "annulus_inner_factor": 0.8}),
                                  (CarPlot, {"lon_min": 0, "lon_max": 10, "lat_min": -1, "lat_max": 1, "d_min": 2, "d_max": 1})]:
            with self.subTest(constructor=constructor.__name__, data=data), self.assertRaises(ValidationError):
                constructor(**data)
        for csv in ["l,b,d\nhello,0,1\n", "l,b,d\n1,0,inf\n"]:
            _, error = science.validate_and_prepare_dataframe(pd.read_csv(io.StringIO(csv)))
            self.assertIsNotNone(error)

    def test_schematic_is_a_png_artifact(self):
        with tempfile.TemporaryDirectory() as directory:
            settings = SimpleNamespace(output_dir=Path(directory), dustmaps_fits_path="", public_base_url="")
            result = execute("schematic", BubbleSchematic(diameter=1).model_dump(), settings)
            self.assertTrue((Path(directory) / result["filename"]).read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
            self.assertEqual(result["url"], "/files/" + result["filename"])


if __name__ == "__main__":
    unittest.main()
