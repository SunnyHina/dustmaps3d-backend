import tempfile
import unittest
from pathlib import Path

from pydantic import ValidationError

from dustmaps_backend.toolkit import ExtinctionInput, calculate_extinction
from dustmaps_backend.visibility import VisibilityInput, calculate_visibility, utc_interval


class AuxiliaryTests(unittest.TestCase):
    def test_visibility_validation_timezone_and_real_png(self):
        fields = dict(start_date="2026-10-08", end_date="2026-10-08",
                      observatory="xinglong", target="Polaris")
        params = VisibilityInput(**fields)
        start, end = utc_interval(params)
        self.assertEqual(start.isoformat(), "2026-10-07T16:00:00+00:00")
        self.assertEqual((end - start).total_seconds(), 86400)
        dst = VisibilityInput(**{**fields, "start_date": "2026-03-08", "end_date": "2026-03-08",
                                 "timezone": "America/New_York"})
        begin, finish = utc_interval(dst)
        self.assertEqual((finish - begin).total_seconds(), 23 * 3600)
        for invalid in ({"end_date": "2026-10-15"}, {"end_date": "2026-10-07"},
                        {"timezone": "Invalid/Zone"}, {"target": None, "lon": 400, "lat": 0}):
            with self.assertRaises(ValidationError):
                VisibilityInput(**{**fields, **invalid})
        with tempfile.TemporaryDirectory() as output:
            result = calculate_visibility(params, output)
            self.assertTrue((Path(output) / result["filename"]).read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
            # Polaris remains near this observatory's latitude all day.
            self.assertTrue(all(38 < altitude < 43 for altitude in result["altitude_deg"]))
            self.assertEqual(len(result["altitude_deg"]), 200)

    def test_empirical_coefficient_and_filter_validation(self):
        result = calculate_extinction(ExtinctionInput(use_2023=True, band="Ks"))
        self.assertAlmostEqual(result["result"], 0.306)
        self.assertEqual(result["unit"], "dimensionless")
        with self.assertRaises(ValidationError):
            ExtinctionInput(use_2023=True, band="GAIA3.Gbp")


if __name__ == "__main__":
    unittest.main()
