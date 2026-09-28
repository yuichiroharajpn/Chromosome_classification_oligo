"""Fail on reversed/nonmatching protein IDs or malformed Larger/Smaller annotations."""
import csv
import sys


def rows(path):
    with open(path, newline='') as stream:
        return list(csv.reader(stream, delimiter='\t'))


def validate(pairs, target_bed, reference_bed, *class_files):
    pair_rows = rows(pairs)
    if not pair_rows or any(len(row) != 2 for row in pair_rows):
        raise ValueError('One2OnePairs must contain two tab-separated columns without a header')
    beds = [rows(target_bed), rows(reference_bed)]
    for index, bed in enumerate(beds):
        if not bed or any(len(row) != 4 for row in bed):
            raise ValueError('Canonical protein BEDs must contain exactly four columns')
        identifiers = {row[3] for row in bed}
        if not any(row[index] in identifiers for row in pair_rows):
            raise ValueError('Pair orientation or IDs do not match target/reference BEDs')
    for path in class_files:
        entries = rows(path)
        if not entries or any(len(row) != 2 or row[1] not in ('Larger', 'Smaller') for row in entries):
            raise ValueError(f'{path}: expected two columns, chromosome and Larger/Smaller, no header')
        if len({row[0] for row in entries}) != len(entries):
            raise ValueError(f'{path}: duplicate chromosome annotations')


if __name__ == '__main__':
    if len(sys.argv) != 6:
        sys.exit('Usage: validate_inputs.py PAIRS TARGET_BED REFERENCE_BED REFERENCE_CLASSES ANCESTRAL_CLASSES')
    try:
        validate(*sys.argv[1:])
    except (OSError, ValueError) as error:
        sys.exit(str(error))
