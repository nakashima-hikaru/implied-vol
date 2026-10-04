"""Checks for algorithm-change checksums and complete paired benchmark inputs."""

import unittest

import compare_experimental as compare


def csv_output(version="baseline", solver="experimental", rounds=2,
               experimental_paths=True, as1_paths=False, dpoly_paths=False,
               omit=None, replace=None, extra=None):
    lines = [",".join(compare.FIELDS)]
    for number in range(rounds):
        for case in compare.cases_for(solver, experimental_paths, as1_paths, dpoly_paths):
            if (case, number) == omit:
                continue
            row = [case, str(number), "100.0", "1" if version == "baseline" else "1.1"]
            if (case, number) == ("atm", 0) and replace is not None:
                row = replace
            lines.append(",".join(row))
    if extra is not None:
        lines.append(extra)
    return "\n".join(lines) + "\n"


def all_samples(blocks=2, rounds=2, **flags):
    samples = []
    for block in range(blocks):
        for solver in compare.SOLVERS:
            for position, version in enumerate(compare.version_order(block)):
                samples.extend(compare.parse_samples(csv_output(version, solver, rounds, **flags),
                    version, solver, rounds, block, position, **flags))
    return samples


class CsvTests(unittest.TestCase):
    def test_population_flags_preserve_default_and_require_opt_in_rows(self):
        self.assertEqual(len(compare.cases_for("experimental")), 22)
        self.assertEqual(len(compare.cases_for("hybrid", as1_paths=True, dpoly_paths=True)), 12)
        flags = {"experimental_paths": False, "as1_paths": True, "dpoly_paths": True}
        samples = compare.parse_samples(csv_output(**flags), "baseline", "experimental", 2, **flags)
        self.assertEqual(len(samples), 2 * 19)
        with self.assertRaisesRegex(compare.ComparisonError, "missing case/round"):
            compare.parse_samples(csv_output(omit=("as1_mixed", 1), **flags),
                                  "baseline", "experimental", 2, **flags)

    def test_missing_duplicate_and_unknown_case_or_round_fail(self):
        for options, error in [
            ({"omit": ("mixed_full", 1)}, "missing case/round"),
            ({"extra": "atm,0,100,1"}, "duplicate round"),
            ({"replace": ["unknown", "0", "100", "1"]}, "unknown case"),
            ({"replace": ["atm", "2", "100", "1"]}, "unexpected round"),
        ]:
            with self.subTest(options=options), self.assertRaisesRegex(compare.ComparisonError, error):
                compare.parse_samples(csv_output(**options), "baseline", "experimental", 2)

    def test_invalid_timing_and_checksum_fail(self):
        for value in ["nan", "inf", "-inf", "0", "-1", "bad"]:
            for row in [["atm", "0", value, "1"], ["atm", "0", "100", value]]:
                with self.subTest(row=row), self.assertRaises(compare.ComparisonError):
                    compare.parse_samples(csv_output(replace=row), "baseline", "experimental", 2)

    def test_invalid_header_and_column_count_fail(self):
        for output in ["", csv_output().replace("ns_per_call", "ns"),
                       csv_output(extra="atm,0,100,1,extra"), csv_output(extra="")]:
            with self.subTest(output=output), self.assertRaises(compare.ComparisonError):
                compare.parse_samples(output, "baseline", "experimental", 2)


class PairedAggregationTests(unittest.TestCase):
    def test_cross_version_changed_checksum_is_allowed_and_recorded(self):
        summary = compare.summarize(all_samples(), 2, 2)
        self.assertEqual(len(summary), 22 + 12)
        self.assertTrue(all(row["samples_per_version"] == 8 for row in summary))
        self.assertTrue(all(row["checksum_changed_between_versions"] for row in summary))
        self.assertTrue(all(row["baseline_checksum"] == "1" and row["candidate_checksum"] == "1.1"
                            for row in summary))

    def test_within_version_checksum_drift_is_rejected(self):
        samples = all_samples()
        samples[0]["checksum"] = "1.01"
        with self.assertRaisesRegex(compare.ComparisonError, "checksum changed within version"):
            compare.summarize(samples, 2, 2)

    def test_missing_and_duplicate_invocation_samples_fail(self):
        samples = all_samples()
        with self.assertRaisesRegex(compare.ComparisonError, "missing samples"):
            compare.summarize(samples[1:], 2, 2)
        with self.assertRaisesRegex(compare.ComparisonError, "duplicate invocation/round"):
            compare.summarize(samples + [samples[0]], 2, 2)

    def test_wrong_abba_position_fails(self):
        samples = all_samples()
        samples[0]["position"] = 1
        with self.assertRaisesRegex(compare.ComparisonError, "position violates ABBA"):
            compare.summarize(samples, 2, 2)

    def test_block_comparisons_and_median_change_use_candidate_over_baseline(self):
        samples = all_samples()
        for row in samples:
            if row["version"] == "candidate":
                row["ns"] = 75.0
        summary = compare.summarize(samples, 2, 2)
        self.assertTrue(all(row["change_percent"] == -25.0 for row in summary))
        self.assertTrue(all(row["paired_block_changes_percent"] == [-25.0, -25.0] for row in summary))


if __name__ == "__main__":
    unittest.main()
