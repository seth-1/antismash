# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from io import StringIO
import unittest

from Bio import SearchIO

from antismash.common import path
from antismash.modules.active_site_finder import common

def build_hmmpfam_results():
    return SearchIO.parse(open(path.get_full_path(__file__, 'data', 'hmmpfam.results')), "hmmer2-text")

class TestCommon(unittest.TestCase):
    def test_get_signature(self):
        query = "TYLVTGGAGGIGGQLALWLAD-QGARHLLLTGRS-A-L--PEQdavvsethpqaTAVAVLRQLRERGVNVTYKAVDVADAHAMQATLESRRRA-GM--PPVRGVFHAAGVIDYTLLSDMSGAEMDRVLAAKVSGAWNLHRLLR-EES----VEAFVLFSSGSALLSSPMLGGYAAGNAFLDALAHHRHAQGL--SGTVVNWGFWD--"
        aln = "tYLitGGlGGLGlslArWLaerrGARrLvLlSRslglpllpsp...........eaqellaeLealGarVrvvacDVtdraavrrllaeiraldtlespPirGViHaAgVLrDallenmtaedfrrVlaPKVdGawnLHeatreddppegsLDFFvlFSSiagllGnpGQanYAAANaFLDAlAryRRarGLRGpAlsinWGaWadv"
        positions = [102]

        assert common.get_signature(query, aln, positions) == "Y"
        assert common.get_signature(query, aln, positions, expected=None) == "Y"

#        queries = list(build_hmmpfam_results())
#        assert len(queries) == 2
#        print(dir(queries[0]))
#        assert queries[0].aln




PK_KR_HMMPFAM_RESULTS = """HMMER 2.3.2 (Oct 2003)
Copyright (C) 1992-2003 HHMI/Washington University School of Medicine
Freely distributed under the GNU General Public License (GPL)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
HMM file:                 antismash/modules/active_site_finder/data/PKSI-KR.hmm2
Sequence file:            -
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query sequence: feature0
Accession:      [none]
Description:    [none]

Scores for sequence family classification (score includes all domains):
Model    Description                                    Score    E-value  N
-------- -----------                                    -----    ------- ---
PKSI-KR                                                 273.3    5.4e-83   1

Parsed for domains:
Model    Domain  seq-f seq-t    hmm-f hmm-t      score  E-value
-------- ------- ----- -----    ----- -----      -----  -------
PKSI-KR    1/1       2   191 .]     1   196 []   273.3  5.4e-83

Alignments of top-scoring domains:
PKSI-KR: domain 1 of 1, from 2 to 191: score 273.3, E = 5.4e-83
                   *->tYLitGGlGGLGlslArWLaerrGARrLvLlSRslglpllpsp....
                      tYL+tGG+GG+G  lA WLa+ +GAR+L L++Rs +    p ++
    feature0     2    TYLVTGGAGGIGGQLALWLAD-QGARHLLLTGRS-A-L--PEQdavv 43

                   .......eaqellaeLealGarVrvvacDVtdraavrrllaeiraldtle
                   ++++++  a + l++L+++G++V++ a+DV+d++a++++l++ r+  ++
    feature0    44 sethpqaTAVAVLRQLRERGVNVTYKAVDVADAHAMQATLESRRRA-GM- 91

                   spPirGViHaAgVLrDallenmtaedfrrVlaPKVdGawnLHeatreddp
                    pP+rGV+HaAgV++  ll++m+ ++++rVla+KV GawnLH+++r +++
    feature0    92 -PPVRGVFHAAGVIDYTLLSDMSGAEMDRVLAAKVSGAWNLHRLLR-EES 139

                   pegsLDFFvlFSSiagllGnpGQanYAAANaFLDAlAryRRarGLRGpAl
                       ++ FvlFSS ++ll +p  + YAA+NaFLDAlA++R+a+GL   ++
    feature0   140 ----VEAFVLFSSGSALLSSPMLGGYAAGNAFLDALAHHRHAQGL--SGT 183

                   sinWGaWadv<-*
                   ++nWG W
    feature0   184 VVNWGFWD--    191

//"""
