"""Contract checks for the volfi same-executable A/B benchmark."""
import unittest
import compare_volfi as bench


class ParseTests(unittest.TestCase):
    def test_parse_and_missing_rows(self):
        header = ",".join(bench.FIELDS)
        output = header + "\nbook_a,0,100.0,5.0,1,1\nbook_a,1,120.0,5.0,1,1\n"
        rows = bench.parse_output(output, "book_a", "experimental", 2)
        self.assertEqual([x["round"] for x in rows], [0, 1])
        with self.assertRaises(ValueError):
            bench.parse_output(output[:-25], "book_a", "experimental", 2)

    def test_bad_rows(self):
        header = ",".join(bench.FIELDS) + "\n"
        for row in ("book_b,0,100,5,1,1", "book_a,0,nan,5,1,1",
                    "book_a,0,100,5,2,1", "book_a,0,100,5,1,1,EXTRA"):
            with self.subTest(row=row), self.assertRaises(ValueError):
                bench.parse_output(header + row + "\n", "book_a", "experimental", 1)

    def test_rotation(self):
        self.assertEqual(bench.rotation(0), ("experimental", "volfi-book", "volfi-book", "experimental"))
        self.assertEqual(bench.rotation(1), ("volfi-book", "experimental", "experimental", "volfi-book"))


if __name__ == "__main__":
    unittest.main()
