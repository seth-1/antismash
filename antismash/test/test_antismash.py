# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import unittest

from antismash import main
from antismash.config.args import build_parser

class TestAntismash(unittest.TestCase):
    def setUp(self):
        self.all_modules = main.get_all_modules()
        self.default_options = build_parser(modules=self.all_modules).parse_args([])

    def test_default_options(self):
        # default options should work with all normal modules
        options = self.default_options
        assert main.verify_options(options, self.all_modules)

        # adding an incompatibility should not be ok
        options.tta = True
        options.input_type = 'prot'
        assert not main.verify_options(options, self.all_modules)

        # and just to be sure the above isn't just because tta
        options.input_type = 'nucl'
        assert main.verify_options(options, self.all_modules)

    def test_help_options(self):
        for option in ["--list-plugins", "--check-prereqs"]:
            options = build_parser(modules=self.all_modules).parse_args([option])
            ret_val = main.run_antismash("", options)
            assert ret_val == 0
