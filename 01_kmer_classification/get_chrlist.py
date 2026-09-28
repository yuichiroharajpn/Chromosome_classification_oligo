#!/usr/bin/env python3
"""Translate the original assembly-report chromosome selection literally."""
import re
import sys


def chromosome_lines(report):
    accession_column = 6 if "GCF" in str(report) else 4
    with open(report, encoding="utf-8", newline="") as source:
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            name, role, molecule = fields[:3]
            if role != "assembled-molecule" or "MT" in name or "MT" in molecule:
                continue
            if molecule.startswith("LG"):
                chromosome = molecule
            elif (re.search(r"^[0-9XYZWIV]+[ABCDLSpqab]?$", molecule)
                  or re.search(r"8_10|9_10[LS]|19_9|3_X|16_21|1A[ab]|^LGE22|^LG25[ab]|^[A-F]\d+|^[a-h]$|25LG[12]|\.part[012]$|^22\.[12]$", molecule)):
                chromosome = ("LG" if name.startswith("LG") else "chr") + molecule
            else:
                chromosome = molecule.replace("align_Mm", "chr", 1)
            accession = fields[accession_column]
            if accession != "na":
                yield f"{chromosome}\t{accession}\t{fields[8]}\n"


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: get_chrlist.py ASSEMBLY_REPORT.txt")
    sys.stdout.writelines(chromosome_lines(sys.argv[1]))
