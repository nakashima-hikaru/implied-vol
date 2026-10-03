"""Contract checks for Hybrid/Jaeckel/C++ benchmark CSV and report aggregation."""

import unittest

import compare_solvers as compare


def csv_output(rounds=2, solver="hybrid", replace=None, omit=None, extra=None):
    lines = [",".join(compare.FIELDS)]
    for round_number in range(rounds):
        for case in compare.CASES:
            key = (case, round_number)
            if key == omit:
                continue
            row = [case, str(round_number), "100.0", str(1 + compare.SOLVERS.index(solver))]
            if key == ("atm", 0) and replace is not None:
                row = replace
            lines.append(",".join(row))
    if extra is not None:
        lines.append(extra)
    return "\n".join(lines) + "\n"


def all_samples(blocks=2, rounds=2):
    return [row for block in range(blocks) for solver in compare.SOLVERS
            for row in compare.parse_samples(csv_output(rounds, solver), solver, rounds, block)]


class CsvContractTests(unittest.TestCase):
    def test_complete_csv_preserves_all_rounds_and_cases(self):
        rows = compare.parse_samples(csv_output(), "hybrid", 2)
        self.assertEqual(len(rows), 2 * len(compare.CASES))
        self.assertEqual({row["case"] for row in rows}, set(compare.CASES))

    def test_missing_case_round_fails(self):
        with self.assertRaisesRegex(compare.ComparisonError, "missing case/round"):
            compare.parse_samples(csv_output(omit=("mixed_full", 1)), "hybrid", 2)

    def test_duplicate_round_fails(self):
        with self.assertRaisesRegex(compare.ComparisonError, "duplicate round"):
            compare.parse_samples(csv_output(extra="atm,0,100,1"), "hybrid", 2)

    def test_unexpected_round_fails(self):
        with self.assertRaisesRegex(compare.ComparisonError, "unexpected round"):
            compare.parse_samples(csv_output(replace=["atm", "2", "100", "1"]), "hybrid", 2)

    def test_unknown_case_fails(self):
        with self.assertRaisesRegex(compare.ComparisonError, "unknown case"):
            compare.parse_samples(csv_output(replace=["extra", "0", "100", "1"]), "hybrid", 2)

    def test_nonfinite_or_nonpositive_timing_fails(self):
        for timing in ["nan", "inf", "-inf", "0", "-1"]:
            with self.subTest(timing=timing), self.assertRaises(compare.ComparisonError):
                compare.parse_samples(csv_output(replace=["atm", "0", timing, "1"]), "hybrid", 2)

    def test_nonfinite_checksum_fails(self):
        for checksum in ["nan", "inf", "-inf"]:
            with self.subTest(checksum=checksum), self.assertRaises(compare.ComparisonError):
                compare.parse_samples(csv_output(replace=["atm", "0", "100", checksum]), "hybrid", 2)

    def test_bad_header_or_column_count_fails(self):
        for output in ["", csv_output().replace("ns_per_call", "ns"),
                       csv_output(extra="atm,0,100,1,extra"), csv_output(extra="")]:
            with self.subTest(output=output), self.assertRaises(compare.ComparisonError):
                compare.parse_samples(output, "hybrid", 2)


class AggregationTests(unittest.TestCase):
    def test_cross_solver_checksum_differences_are_allowed(self):
        summary = compare.summarize(all_samples(), 2, 2)
        self.assertEqual(len(summary), len(compare.SOLVERS) * len(compare.CASES))
        self.assertTrue(all(row["samples"] == 4 for row in summary))

    def test_checksum_drift_within_solver_fails(self):
        samples = all_samples()
        samples[0]["checksum"] = "1.1"
        with self.assertRaisesRegex(compare.ComparisonError, "checksum changed within solver"):
            compare.summarize(samples, 2, 2)

    def test_missing_solver_case_sample_fails(self):
        with self.assertRaisesRegex(compare.ComparisonError, "missing samples"):
            compare.summarize(all_samples()[1:], 2, 2)

    def test_duplicate_block_round_fails(self):
        samples = all_samples()
        with self.assertRaisesRegex(compare.ComparisonError, "duplicate block/round"):
            compare.summarize(samples + [samples[0]], 2, 2)

    def test_relative_percent_uses_cpp_time(self):
        samples = all_samples()
        for row in samples:
            if row["solver"] == "hybrid":
                row["ns_per_call"] = 75.0
        summary = compare.summarize(samples, 2, 2)
        percentages = {row["time_change_vs_cpp_percent"] for row in summary
                       if row["solver"] == "hybrid"}
        self.assertEqual(percentages, {-25.0})

    def test_order_rotates_and_reverses_without_omitting_solvers(self):
        self.assertEqual(compare.measurement_order(0), compare.SOLVERS)
        self.assertEqual(compare.measurement_order(1),
                         tuple(reversed(compare.SOLVERS[1:] + compare.SOLVERS[:1])))
        for block in range(10):
            self.assertCountEqual(compare.measurement_order(block), compare.SOLVERS)


if __name__ == "__main__":
    unittest.main()
