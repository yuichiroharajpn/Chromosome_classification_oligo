"""Exercise shell orchestration without installing the scientific tools."""
import gzip
import json
import os
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

WRAPPER = Path(__file__).resolve().parents[1] / 'run_01_kmer_classification.sh'
MOCK = '''#!PYTHON
import json, os, pathlib, sys
name = pathlib.Path(sys.argv[0]).name
args = sys.argv[1:]
with open(os.environ['CALL_LOG'], 'a') as log:
    log.write(json.dumps([name, args]) + '\\n')
if name == 'bgzip':
    if os.environ.get('FAIL_BGZIP'): sys.exit(7)
    sys.stdout.buffer.write(sys.stdin.buffer.read())
elif name == 'samtools' and len(args) == 3:
    print('>' + args[2] + '\\nACGTACGT')
elif name == 'jellyfish':
    if args[0] == 'count':
        sys.stdin.buffer.read()
        if os.environ.get('FAIL_COUNT') and args[args.index('-m') + 1] == '4': sys.exit(9)
        pathlib.Path(args[args.index('-o') + 1]).write_bytes(b'database')
    else: print('AC\\t2')
'''.replace('PYTHON', sys.executable)


class Wrapper(unittest.TestCase):
    def execute(self, fail=False, mode=None, cores=1, fail_count=False):
        with tempfile.TemporaryDirectory(prefix='kmer test ') as directory:
            root = Path(directory)
            bin_dir = root / 'bin'
            bin_dir.mkdir()
            for name in ('bgzip', 'samtools', 'jellyfish', 'Rscript'):
                tool = bin_dir / name
                tool.write_text(MOCK)
                tool.chmod(0o755)
            if os.environ.get('PARALLEL_TEST_EXEC'):
                (bin_dir / 'parallel').symlink_to(os.environ['PARALLEL_TEST_EXEC'])
            assembly = root / 'GCF_example_genomic.fna.gz'
            with gzip.open(assembly, 'wb') as stream:
                stream.write(b'>acc\nACGTACGT\n')
            assembly.with_name('GCF_example_assembly_report.txt').write_text(
                'seq\tassembled-molecule\t1\tChromosome\tgb\t=\tacc\tPrimary Assembly\t101\tna\n')
            log = root / 'calls.jsonl'
            env = dict(os.environ, PATH=str(bin_dir) + os.pathsep + os.environ['PATH'], CALL_LOG=str(log))
            if fail:
                env['FAIL_BGZIP'] = '1'
            if fail_count:
                env['FAIL_COUNT'] = '1'
            output = root / 'results'
            result = subprocess.run(['bash', str(WRAPPER), str(assembly), '--output-dir', str(output), '--cores', str(cores)] + ([mode] if mode else []), env=env, capture_output=True)
            calls = [json.loads(line) for line in log.read_text().splitlines()]
            if fail_count:
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(any(name == 'Rscript' and args[0] != '-e' for name, args in calls))
                return
            if fail:
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(any(name == 'samtools' for name, args in calls))
                return
            self.assertEqual(result.returncode, 0, result.stderr.decode())
            self.assertEqual((output / 'chromosome_list.txt').read_bytes(), b'chr1\tacc\t101\n')
            expected_regions = ['acc'] if mode is None else (['acc', 'acc:5-95'] if mode == '--create-mid' else ['acc:5-95'])
            expected_files = ['chr1.fa.gz'] if mode is None else (['chr1.fa.gz', 'chr1_mid.fa.gz'] if mode == '--create-mid' else ['chr1_mid.fa.gz'])
            self.assertEqual(sorted(p.name for p in (output / 'chromosomes').glob('*.fa.gz')), expected_files)
            self.assertEqual([args[2] for name, args in calls if name == 'samtools' and len(args) == 3], expected_regions)
            counts = [args for name, args in calls if name == 'jellyfish' and args[0] == 'count']
            self.assertEqual(sorted(int(args[2]) for args in counts), sorted(list(range(2, 10)) * len(expected_regions)))
            self.assertEqual(len(list((output / 'chromosomes').glob('*mer.txt'))), 8 * len(expected_regions))
            analyses = [Path(args[0]).name for name, args in calls if name == 'Rscript' and args[0] != '-e']
            self.assertEqual(analyses, ['summary_kmerfreq.R'] * 8 + ['analyze_freq.R'] * 8 + ['oligo_coordinates.R', 'top_oligos.R'])
            r_calls = [args for name, args in calls if name == 'Rscript' and args[0] != '-e']
            prefix, suffix = ('', '') if mode is None else ('mid_', '_mid')
            self.assertEqual(sorted(args[1:] for args in r_calls[:8]), [[str(k), f'{prefix}{k}mer', f'{k}mer{suffix}'] for k in range(2, 10)])
            self.assertEqual(sorted(args[1] for args in r_calls[8:16]), [f'{k}mer{suffix}_relfreq.txt' for k in range(2, 10)])
            self.assertEqual(r_calls[-2][1], f'6mer{suffix}_relfreq.txt')
            self.assertEqual(r_calls[-1][1], f'6mer{suffix}_relfreq.coord.txt')
            self.assertTrue(assembly.exists())

    def test_invalid_cores(self):
        for value in ('0', '9', '-1', '1.5', 'abc', '08'):
            result = subprocess.run(['bash', str(WRAPPER), '--cores', value], capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(b'--cores must be an integer from 1 to 8', result.stderr)

    @unittest.skipUnless(os.environ.get('PARALLEL_TEST_EXEC') or shutil.which('parallel'), 'GNU parallel unavailable')
    def test_parallel_modes(self):
        for cores, mode in ((2, None), (4, '--create-mid'), (8, '--create-mid-only')):
            with self.subTest(cores=cores, mode=mode):
                self.execute(mode=mode, cores=cores)

    @unittest.skipUnless(os.environ.get('PARALLEL_TEST_EXEC') or shutil.which('parallel'), 'GNU parallel unavailable')
    def test_parallel_count_failure(self):
        self.execute(cores=2, fail_count=True)

    def test_complete_orchestration(self):
        self.execute()

    def test_create_mid(self):
        self.execute(mode='--create-mid')

    def test_create_mid_only(self):
        self.execute(mode='--create-mid-only')

    def test_conflicting_options(self):
        result = subprocess.run(['bash', str(WRAPPER), '--create-mid', '--create-mid-only'], capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b'mutually exclusive', result.stderr)

    def test_upstream_failure_stops_run(self):
        self.execute(fail=True)


if __name__ == '__main__':
    unittest.main()
