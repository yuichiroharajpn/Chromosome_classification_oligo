"""Regression checks against the private legacy source, when available."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

PUBLIC = Path(__file__).resolve().parents[1]
LEGACY = PUBLIC.parents[1] / "tmp" / "01_kmer_classification"
sys.path.insert(0, str(PUBLIC))


class Compatibility(unittest.TestCase):
    @unittest.skipUnless(LEGACY.exists(), "legacy source not distributed")
    def test_report_bytes(self):
        molecules = ["1", "X", "W", "IV", "1A", "LG2", "8_10", "9_10S", "19_9", "3_X", "16_21", "1Ab", "LGE22", "LG25a", "A12", "a", "25LG1", "x.part0", "22.1", "align_Mm2", "foo", "MT", "1MT", "01", "1q"]
        with tempfile.TemporaryDirectory() as temporary:
            for prefix in ("GCF", "GCA"):
                report = Path(temporary) / (prefix + "_assembly_report.txt")
                records = ["# comment\n"]
                for index, molecule in enumerate(molecules):
                    for name in ("seq", "LGseq", "MTseq"):
                        for role in ("assembled-molecule", "unlocalized-scaffold"):
                            for accession in ("RS1", "na"):
                                records.append("\t".join([name, role, molecule, "Chromosome", accession, "=", accession, "Primary Assembly", str(1001 + index), "na"]) + "\n")
                report.write_text("".join(records))
                expected = subprocess.check_output(["perl", str(LEGACY / "get_chrlist.pl"), str(report)])
                actual = subprocess.check_output([sys.executable, str(PUBLIC / "get_chrlist.py"), str(report)])
                self.assertEqual(actual, expected)

    def test_midpoint_rounding(self):
        lengths = list(range(1, 1001)) + [100000001, 2147483647, 4294967295]
        expected = subprocess.check_output(["perl", "-e", 'for (@ARGV) {printf "%d-%d\\n", $_*0.05, $_*0.95}', *map(str, lengths)]).decode().splitlines()
        actual = subprocess.check_output([
            "bash", "-c", 'source "$1"; shift; for n in "$@"; do mid_region acc "$n"; done',
            "test", str(PUBLIC / "run_01_kmer_classification.sh"), *map(str, lengths)
        ]).decode().replace("acc:", "").splitlines()
        self.assertEqual(actual, expected)

    @unittest.skipUnless(LEGACY.exists(), "legacy source not distributed")
    def test_r_sources_unchanged(self):
        for old, new in [("analysis_freq.v4.R", "analyze_freq.R"), ("analysis_freq.v3.1.R", "oligo_coordinates.R"), ("summary_kmerfreq.R", "summary_kmerfreq.R"), ("getoligo_Dim1top100.R", "top_oligos.R")]:
            self.assertEqual((LEGACY / old).read_bytes(), (PUBLIC / new).read_bytes())


if __name__ == "__main__":
    unittest.main()
