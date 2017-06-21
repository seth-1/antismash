# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import os
import sys
import unittest

from antismash.main import run_antismash, gather_modules
from antismash.config.args import build_parser, Config

class TestAntismash(unittest.TestCase):
    def setUp(self):
        self.old_argv = sys.argv
        sys.argv = ["run_antismash.py"]
        self.parser = build_parser(modules=gather_modules(with_genefinding=True))
        self.default_options = self.parser.parse_args(sys.argv)
        Config(self.default_options)

    def tearDown(self):
        Config().__dict__.clear()
        sys.argv = self.old_argv

    def test_nisin_minimal(self):
        path = os.path.join(os.path.dirname(__file__), "data", "nisin.gbk")
        run_antismash(path, Config({"minimal": True}))
