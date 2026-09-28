"""Compare public table bytes with legacy outputs, allowing only requested changes."""
from pathlib import Path
import subprocess
import os
import sys
import tempfile
import unittest

PUBLIC = Path(__file__).resolve().parents[1]
LEGACY = PUBLIC.parents[1] / 'tmp/02_cross_species_assignment'
sys.path.insert(0, str(PUBLIC))
from select_best import select_best
from validate_inputs import validate


def execute(script, *args, cwd=None):
    return subprocess.check_output([sys.executable, str(script), *map(str, args)], cwd=cwd, stderr=subprocess.PIPE)


def rename_classes(data):
    # Replace complete tabular values only; leave headers and all other bytes intact.
    data = data.replace(b'\tMac\t', b'\tLarger\t').replace(b'\tmin\t', b'\tSmaller\t')
    return data.replace(b'\tMac\t', b'\tLarger\t').replace(b'\tmin\t', b'\tSmaller\t').replace(b'\tMac\n', b'\tLarger\n').replace(b'\tmin\n', b'\tSmaller\n')


class Selection(unittest.TestCase):
    def row(self, k, height, score, kind='AutoXZ'):
        return dict(Kmer=k, height_class=height, Informedness=score, chromosome_class=kind)

    def test_constraints_and_ties(self):
        rows = [self.row(2, 'L1', 1), self.row(3, 'L1', 1, 'Auto'),
                self.row(9, 'L2', .8), self.row(4, 'L3', .8), self.row(4, 'L4', .8)]
        self.assertEqual(select_best(rows), rows[3])

    def test_unrounded_score(self):
        rows = [self.row(3, 'L1', .80001), self.row(4, 'L2', .80002)]
        self.assertEqual(select_best(rows), rows[1])

    def test_no_candidates(self):
        with self.assertRaises(ValueError):
            select_best([self.row(2, 'L1', 1), self.row(3, 'L2', float('nan'))])


class Config(unittest.TestCase):
    def test_invalid_entries(self):
        for content, message in [('unknown = value\n', b'unknown key'),
                                 ('pairs = a\npairs = b\n', b'duplicate key'),
                                 ('pairs value\n', b'expected key = path')]:
            with tempfile.TemporaryDirectory() as directory:
                config=Path(directory)/'config.cfg'
                config.write_text(content)
                result=subprocess.run(['bash',str(PUBLIC/'run_02_cross_species_assignment.sh'),
                                       '--config',str(config)], capture_output=True)
                self.assertNotEqual(result.returncode,0)
                self.assertIn(message,result.stderr)

    def test_shell_text_is_not_executed(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            marker=root/'executed'
            config=root/'config.cfg'
            config.write_text('pairs = $(touch '+str(marker)+')\n')
            subprocess.run(['bash',str(PUBLIC/'run_02_cross_species_assignment.sh'),
                            '--config',str(config)],capture_output=True)
            self.assertFalse(marker.exists())


@unittest.skipUnless(LEGACY.exists(), 'private legacy scripts unavailable')
class Compatibility(unittest.TestCase):
    def test_supplied_class_tables(self):
        data = PUBLIC
        for name in ('chromosome_list.class.txt', 'Huang_ancestral.class.txt'):
            original = (LEGACY/'ref_chick'/name).read_bytes()
            public = (data/'ref_chick'/name).read_bytes()
            self.assertEqual(public, rename_classes(original))
            self.assertEqual({line.split(b'\t')[1] for line in public.splitlines()}, {b'Larger', b'Smaller'})
        validate(data/'inparanoid_out/One2OnePairs',
                 data/'ref_zebrafinch/GCF_003957565.2_bTaeGut1.4.pri_protein.canonical.bed',
                 data/'ref_chick/chicken.v23.pep.canonical.bed',
                 data/'ref_chick/chromosome_list.class.txt',
                 data/'ref_chick/Huang_ancestral.class.txt')

    def test_one_to_one(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'pairs'
            path.write_text('A_ABCDEF\tB_abcdef\nC\tD\nC\tE\nF\tG\nH\tG\nI\tJ\n')
            self.assertEqual(execute(LEGACY/'get1to1.py', path), execute(PUBLIC/'get_one_to_one.py', path))

    def test_r_operations(self):
        # Ignoring comments and whitespace, preserve every R operation in order.
        for old, new in [('macrosynt.R', 'macrosynteny.R'), ('macrosynt.anc.rev.R', 'macrosynteny_ancestral.R')]:
            def statements(path):
                return [line.strip() for line in path.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
            self.assertEqual(statements(LEGACY/old), statements(PUBLIC/new))

    def test_tables(self):
        with tempfile.TemporaryDirectory(prefix='assignment ') as directory:
            root = Path(directory).resolve()
            cluster = 'Chr\tL1\tL2\tL3\tL4\tL5\tLength\nchr1\t1\t1\t1\t1\t1\t100\nchr2\t1\t2\t2\t2\t2\t90\nchr3\t2\t2\t3\t3\t3\t80\nchrZ\t1\t1\t1\t4\t4\t70\nchrW\t2\t2\t2\t5\t5\t60\n'
            for k in (10, 2, 3, 4):
                (root/f'{k}mer_mid_relfreq.dendro.grp.h.txt').write_text(cluster)
            reference = root/'reference.tsv'
            reference.write_text('sp1.Chr\tsp2.Chr\torthologs\tsp1.Chr.ratio\nchr1\tr1\t20\t0.9\nchr2\tr2\t10\t0.8\nchr3\tr3\t9\t0.7\nchrZ\tr1\t8\t0.8\nchrW\tr3\t7\t0.6\n')
            ancestral=root/'ancestral.tsv'
            ancestral.write_text('sp1.Chr\tsp2.Chr\torthologs\nanc1\tchr1\t20\nanc2\tchr2\t10\nanc3\tchr3\t9\nanc4\tchrZ\t8\nanc5\tchrW\t7\n')
            refclasses=root/'refclasses.tsv'
            refclasses.write_text('r1\tMac\nr2\tmin\nr3\tmin\n')
            ancclasses=root/'ancclasses.tsv'
            ancclasses.write_text('anc1\tMac\nanc2\tmin\nanc3\tmin\nanc4\tMac\nanc5\tmin\n')
            public_refclasses = root/'public_refclasses.tsv'
            public_ancclasses = root/'public_ancclasses.tsv'
            public_refclasses.write_bytes(rename_classes(refclasses.read_bytes()))
            public_ancclasses.write_bytes(rename_classes(ancclasses.read_bytes()))
            for old, new, synteny, classes in [('rankdist.ref.h.py','rank_reference.py',reference,refclasses),('rankdist.ref.anc_jv.h.rev.py','rank_ancestral.py',ancestral,ancclasses)]:
                original = execute(LEGACY/old, root, synteny, classes).decode().splitlines()
                expected = ['\t'.join(original[0].replace('Kmeans','height_class').replace('order','chromosome_class').replace('Dist_comp','Informedness').split('\t'))]
                for line in original[1:]:
                    values=line.split('\t')
                    if values[4] != 'AutoXZ':
                        continue
                    values[3]=str(int(float(values[3])))
                    expected.append('\t'.join(values))
                selection=root/'selection.txt'
                actual=execute(PUBLIC/new, root, synteny, public_refclasses if classes == refclasses else public_ancclasses, '--selection-file', selection)
                self.assertEqual(actual, ('\n'.join(expected)+'\n').encode())
                self.assertTrue(selection.read_text().startswith('3\t'))
            path=root/'3mer_mid_relfreq.dendro.grp.h.txt'
            self.assertEqual(rename_classes(execute(LEGACY/'match_dendro.py',path,'L2')),execute(PUBLIC/'match_dendro.py',path,'L2'))
            huang=root/'Birds/Gallus_gallus_T2T/ancestral/Huang_ancestral.txt'
            huang.parent.mkdir(parents=True)
            huang.write_text('anc_jv\tChicken.chromosome\nanc1\t1\nanc2\t2\nanc3\t3\nanc4\tZ\nanc5\tW\n')
            self.assertEqual(rename_classes(execute(LEGACY/'match_anc_chrom.py',path,ancestral,ancclasses,'L2',cwd=root)),execute(PUBLIC/'match_anc_chrom.py',path,ancestral,public_ancclasses,'L2',huang))
            pairs=root/'pairs.tsv'
            pairs.write_text('target1\treference1\n')
            target_bed=root/'target.bed'
            target_bed.write_text('chr1\t1\t20\ttarget1\n')
            reference_bed=root/'reference.bed'
            reference_bed.write_text('r1\t1\t20\treference1\n')
            output=root/'wrapper output'
            env=dict(os.environ, PATH=str(Path(sys.executable).parent)+os.pathsep+os.environ['PATH'])
            subprocess.run(['bash',str(PUBLIC/'run_02_cross_species_assignment.sh'),
                            '--kmer-dir',str(root),'--pairs',str(pairs),
                            '--target-bed',str(target_bed),'--reference-bed',str(reference_bed),
                            '--ancestral-table',str(huang),'--reference-classes',str(public_refclasses),
                            '--ancestral-classes',str(public_ancclasses),'--output-dir',str(output),
                            '--reference-synteny',str(reference),'--ancestral-synteny',str(ancestral)],
                           env=env,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            config = root/'example config.cfg'
            configured_output = root/'configured output'
            settings = dict(kmer_dir=root, pairs=pairs, target_bed=target_bed,
                            reference_bed=reference_bed, ancestral_table=huang,
                            reference_classes=public_refclasses, ancestral_classes=public_ancclasses,
                            output=root/'unused config output', reference_synteny=reference,
                            ancestral_synteny=ancestral)
            config.write_text('\n'.join(key+' = '+str(value.relative_to(root)) for key,value in settings.items())+'\n')
            subprocess.run(['bash', str(PUBLIC/'run_02_cross_species_assignment.sh'),
                            '--output-dir', str(configured_output), '--config', str(config)],
                           cwd='/', env=env, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            for expected_file in output.rglob('*'):
                if expected_file.is_file():
                    self.assertEqual(expected_file.read_bytes(),
                        (configured_output/expected_file.relative_to(output)).read_bytes())
            self.assertFalse((root/'unused config output').exists())
            mock_bin = root/'mock bin'
            mock_bin.mkdir()
            mock_r = mock_bin/'Rscript'
            mock_r.write_text('#!' + sys.executable + '\n' +
                "import json, pathlib, sys\n" +
                "args = sys.argv[1:]\n" +
                "with open(" + repr(str(root/'r_calls.jsonl')) + ", 'a') as log: log.write(json.dumps(args) + '\\n')\n" +
                "if args[0] != '-e':\n" +
                "    source = " + repr(str(ancestral)) + " if pathlib.Path(args[0]).name == 'macrosynteny_ancestral.R' else " + repr(str(reference)) + "\n" +
                "    sys.stdout.buffer.write(pathlib.Path(source).read_bytes())\n")
            mock_r.chmod(0o755)
            computed_output = root/'default script output'
            subprocess.run(['bash', str(PUBLIC/'run_02_cross_species_assignment.sh'),
                            '--kmer-dir', str(root), '--pairs', str(pairs),
                            '--target-bed', str(target_bed), '--reference-bed', str(reference_bed),
                            '--ancestral-table', str(huang), '--reference-classes', str(public_refclasses),
                            '--ancestral-classes', str(public_ancclasses), '--output-dir', str(computed_output)],
                           env=dict(env, PATH=str(mock_bin)+os.pathsep+env['PATH']),
                           check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            import json
            calls = [json.loads(line) for line in (root/'r_calls.jsonl').read_text().splitlines()]
            self.assertEqual(calls[-1], [str(PUBLIC/'macrosynteny_ancestral.R'),
                str(computed_output/'inparanoid_out/One2OnePairs.rev'), str(reference_bed), str(target_bed), str(huang)])
            self.assertEqual(calls[-2], [str(PUBLIC/'macrosynteny.R'),
                str(computed_output/'inparanoid_out/One2OnePairs'), str(target_bed), str(reference_bed)])
            for selection_name, match_name, script in [
                ('reference.best.tsv','match_mid_relfreq.dendro.grp.h.best.txt','match_dendro.py'),
                ('ancestral.best.tsv','match_mid_relfreq.dendro.grp.refanc_jv.h.best.refgalGalT2T.txt','match_anc_chrom.py')]:
                k,height=(output/selection_name).read_text().strip().split('\t')
                self.assertEqual(k,'3')
                selected=root/f'{k}mer_mid_relfreq.dendro.grp.h.txt'
                arguments=[selected,height] if script=='match_dendro.py' else [selected,ancestral,public_ancclasses,height,huang]
                self.assertEqual((output/match_name).read_bytes(),execute(PUBLIC/script,*arguments))



if __name__ == '__main__':
    unittest.main()
