import inspect

import httpimport
import os
import shutil
import tempfile
import timeout_decorator
import time

from pathlib import Path
from unittest import TestCase

# import sys
# sys.path.append('../../islandcompare-cli/')
# import islandcompare as cli
with httpimport.github_repo("brinkmanlab", "islandcompare-cli", "islandcompare", branch="master"):
    import islandcompare as cli


HOUR = 60 * 60
path_prefix = Path('./assets')
workflow_outputs = ("Results", "Newick", "Genomic Islands", 'Stitched genomes', 'Prepared data')
pseudomonas_data = list(path_prefix.glob('Pseudomonas/*.gbk'))
enterococcus_data = list(path_prefix.glob('Enterococcus/*.gbk'))
listeria_data = list(path_prefix.glob('Listeria/*.gbk'))
drafts = list(path_prefix.glob('Draft/*.gbk'))
reference = 'NZ_LN870292_1'  # Pseudomonas


class TestBase(TestCase):
    host = os.environ.get('GALAXY_HOST', 'https://galaxy.islandcompare.ca')
    key = os.environ.get('GALAXY_API_KEY', '')

    def setUp(self) -> None:
        super().setUp()
        self.conn = cli.GalaxyInstance(self.host, self.key)

        self.upload_history = cli.get_upload_history(self.conn)
        if len(self.upload_history.get_datasets()):
            # History isn't fresh, delete and recreate
            cli._retryConnection(self.upload_history.delete, purge=True)
            self.upload_history = cli.get_upload_history(self.conn)

        self.workflow = cli.get_workflow(self.conn)
        self.output_path = Path(tempfile.mkdtemp())

    def tearDown(self) -> None:
        super().tearDown()
        shutil.rmtree(self.output_path)
        # delete all histories
        for history in self.conn.histories.list():
            cli._retryConnection(history.delete, purge=True)

    def confirmOutput(self, outputs):
        self.assertIsNotNone(outputs)
        self.assertSetEqual(set(outputs.keys()), set(workflow_outputs))

    def checkErrors(self, errors):
        self.assertEqual(len(errors), 0)

    def run_workflow(self, data, *args, **kwargs):
        label = inspect.currentframe().f_back.f_code.co_name
        ret, err = cli.round_trip(self.upload_history, data, self.workflow, label, self.output_path, *args, **kwargs)
        self.checkErrors(err)
        self.confirmOutput(ret)
        return ret, err


class TestBasic(TestBase):
    @timeout_decorator.timeout(4 * HOUR)
    def test_pair(self):
        self.run_workflow(pseudomonas_data[:2])

    @timeout_decorator.timeout(4 * HOUR)
    def test_pseudomonas(self):
        self.run_workflow(pseudomonas_data)

    @timeout_decorator.timeout(4 * HOUR)
    def test_listeria(self):
        self.run_workflow(listeria_data)

    @timeout_decorator.timeout(30 * HOUR)
    def test_enterococcus_x1023(self):
        self.assertEqual(len(enterococcus_data), 1023)
        self.run_workflow(enterococcus_data)


class TestDraft(TestBase):
    @timeout_decorator.timeout(4 * HOUR)
    def test_draft_noref(self):
        self.run_workflow(drafts)

    @timeout_decorator.timeout(4 * HOUR)
    def test_draft_ref(self):
        self.run_workflow(drafts, reference_id=reference)

    @timeout_decorator.timeout(20 * HOUR)
    def test_mcm_lockup_ERR388703(self):
        self.run_workflow([Path.joinpath(path_prefix, 'crash_mcm', 'ERR388703.gbk'), drafts[0]], reference_id='NZ_LN999987_1')

    @timeout_decorator.timeout(20 * HOUR)
    def test_mcm_lockup_VRE0027(self):
        self.run_workflow([Path.joinpath(path_prefix, 'crash_mcm', 'VRE-0027.gbk'), drafts[0]], reference_id='NZ_LR135364_1')


class TestNewick(TestBase):
    drafts = list(path_prefix.glob('Genomes_With_Trees/Genomes/*.gbk'))
    accession = Path.joinpath(path_prefix, 'Genomes_With_Trees', 'Trees', 'ByAccession.newick')
    filename = Path.joinpath(path_prefix, 'Genomes_With_Trees', 'Trees', 'ByFilename.newick')
    reference = 'NC_010084_1'

    @timeout_decorator.timeout(4 * HOUR)
    def test_by_accession(self):
        self.run_workflow(self.drafts, newick=self.accession, accession=True, reference_id=reference)

    @timeout_decorator.timeout(4 * HOUR)
    def test_by_filename(self):
        self.run_workflow(self.drafts, newick=self.filename, accession=False, reference_id=reference)


class TestBlastCheck(TestBase):
    path_prefix = Path.joinpath(path_prefix, 'Genomes_With_Trees', 'Genomes')

    @timeout_decorator.timeout(4 * HOUR)
    def test_empty(self):
        self.run_workflow([Path.joinpath(self.path_prefix, 'VC14135.gbk'), Path.joinpath(self.path_prefix, 'VC13213.gbk')], reference_id=reference)


class TestCurated(TestBase):
    path_prefix = Path.joinpath(path_prefix, 'Curated')

    @timeout_decorator.timeout(4 * HOUR)
    def test_basic_salmonella(self):
        self.run_workflow([Path.joinpath(self.path_prefix, 'Salmonella', 'CP091999.1.gbk'), Path.joinpath(self.path_prefix, 'Salmonella', 'LT571437.1.gbk')], reference_id=reference)
