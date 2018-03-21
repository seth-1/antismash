# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from collections import defaultdict
import unittest

import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from antismash.common.test import helpers
from antismash.common import secmet
from antismash.common.secmet.feature import Cluster, CDSFeature, Feature, GeneFunction, ClusterBorder, PFAMDomain
from antismash.common.secmet.record import _build_products_from_borders, Record


class TestConversion(unittest.TestCase):
    def test_record_conversion_from_biopython(self):
        before = list(Bio.SeqIO.parse(helpers.get_path_to_nisin_genbank(), "genbank"))[0]
        # sort notes, because direct comparisons otherwise are awful
        for feature in before.features:
            if "note" in feature.qualifiers:
                feature.qualifiers["note"] = sorted(feature.qualifiers["note"])
        before_features = sorted(map(str, before.features))
        type_counts = defaultdict(lambda: 0)
        for feature in before.features:
            type_counts[feature.type] += 1
        record = Record.from_biopython(before, taxon="bacteria")
        after = record.to_biopython()

        # ensure new features are correct
        assert len(before_features) == len(after.features)
        for bef, aft in zip(before_features, sorted(map(str, after.features))):
            assert bef == aft

        # ensure we haven't changed the original record or feature list
        assert id(before) != id(after)
        assert id(before.features) != id(after.features)
        for i in range(len(before.features)):
            assert id(before.features[i]) != id(after.features[i])
        for bef, aft in zip(before_features, sorted(map(str, before.features))):
            assert bef == aft

        # ensure that the counts of each match
        assert type_counts["CDS"] == len(record.get_cds_features())
        assert type_counts["PFAM_domain"] == len(record.get_pfam_domains())
        assert type_counts["cluster"] == len(record.get_clusters())
        assert type_counts["aSDomain"] == len(record.get_antismash_domains())

    def test_pfam_domain(self):
        original = PFAMDomain(FeatureLocation(2, 5), description="test",
                              domain="domainname")
        original.db_xref.append("test-ref")
        original.tool = "toolname"
        original.domain_id = "domain_id"
        original.database = "db"
        original.detection = "someprogram"
        original.evalue = 1e-5
        original.score = 5.
        original.locus_tag = "locus"
        original.label = "somelabel"
        original.translation = "ARNDCQ"
        new = PFAMDomain.from_biopython(original.to_biopython()[0])
        for slot in ["db_xref", "tool", "domain_id", "database", "detection",
                     "evalue", "score", "locus_tag", "label", "translation"]:
            assert getattr(original, slot) == getattr(new, slot)


class TestRecord(unittest.TestCase):
    def test_cluster_numbering(self):
        record = Record(Seq("A"*1000))
        for start, end in [(50, 100), (10, 40), (700, 1000), (0, 9)]:
            cluster = helpers.DummyCluster(start, end)
            record.add_cluster(cluster)
        for i, cluster in enumerate(sorted(list(record.get_clusters()))):
            assert cluster.get_cluster_number() == i + 1

    def test_overlapping_clusters(self):
        record = Record(seq="A"*40)
        record.add_cluster(Cluster(FeatureLocation(10, 40), 0, 0, []))
        with self.assertRaises(ValueError):
            record.add_cluster(Cluster(FeatureLocation(0, 11), 0, 0, []))
        # ok, since ends aren't inclusive
        record.add_cluster(Cluster(FeatureLocation(0, 10), 0, 0, []))

    def test_cds_cluster_linkage(self):
        record = Record("A"*200)
        for start, end in [(50, 100), (10, 90), (0, 9), (150, 200)]:
            record.add_cds_feature(helpers.DummyCDS(start, end))
        for start, end in [(10, 120), (5, 110), (10, 160), (45, 200)]:
            record.clear_clusters()
            cluster = helpers.DummyCluster(start, end)
            record.add_cluster(cluster)
            assert len(cluster.cds_children) == 2
            for cds in cluster.cds_children:
                assert cds.overlaps_with(cluster)

    def test_cds_removal(self):
        record = Record(Seq("A" * 1000))
        cluster = helpers.DummyCluster(0, 1000)
        record.add_cluster(cluster)

        first_cds = helpers.DummyCDS(0, 100, locus_tag="A")
        second_cds = helpers.DummyCDS(200, 300, locus_tag="B")
        record.add_cds_feature(first_cds)
        record.add_cds_feature(second_cds)

        assert len(record.get_cds_features()) == 2
        assert len(cluster.cds_children) == 2

        record.remove_cds_feature(first_cds)

        assert len(record.get_cds_features()) == 1
        assert len(cluster.cds_children) == 1
        assert record.get_cds_features()[0] is list(cluster.cds_children)[0]
        assert record.get_cds_features()[0].locus_tag == "B"

    def test_orphaned_cluster_number(self):
        record = Record(Seq("A" * 1000))
        cluster = helpers.DummyCluster(0, 1000)

        with self.assertRaisesRegex(ValueError, "Cluster not contained in record"):
            print(record.get_cluster_number(cluster))

        with self.assertRaisesRegex(ValueError, "Cluster not contained in record"):
            print(cluster.get_cluster_number())

    def test_products_from_borders(self):
        location = FeatureLocation(1, 10)
        border1 = ClusterBorder(location, "toolA", product="A")
        assert border1.high_priority_product
        border2 = ClusterBorder(location, "toolB", product="B")
        assert border1.high_priority_product
        assert _build_products_from_borders([border1, border2]) == ["A", "B"]
        assert _build_products_from_borders([border2, border1]) == ["B", "A"]
        border2 = ClusterBorder(location, "toolB", product="B", high_priority_product=False)
        assert _build_products_from_borders([border1, border2]) == ["A"]
        assert _build_products_from_borders([border2, border1]) == ["A"]
        border1.high_priority_product = False
        assert _build_products_from_borders([border1, border2]) == ["A", "B"]
        assert _build_products_from_borders([border2, border1]) == ["B", "A"]


class TestFeature(unittest.TestCase):
    def test_overlaps_with(self):
        # no overlap
        feature = helpers.DummyFeature(5, 10)
        assert isinstance(feature, Feature)  # just to be sure it works the way we want
        other = helpers.DummyFeature(100, 110)
        assert not feature.overlaps_with(other) and not other.overlaps_with(feature)
        # completely within
        other = helpers.DummyFeature(0, 20)
        assert feature.overlaps_with(other) and other.overlaps_with(feature)
        # partially within
        other = helpers.DummyFeature(0, 8)
        assert feature.overlaps_with(other) and other.overlaps_with(feature)
        # borders touching
        other = helpers.DummyFeature(0, 5)
        assert not (feature.overlaps_with(other) or other.overlaps_with(feature))

    def test_is_contained_by(self):
        # same location is considered to be contained
        feature = helpers.DummyFeature(5, 10)
        assert feature.is_contained_by(feature)
        for strand in (-1, 1):
            # no overlap
            other = helpers.DummyFeature(15, 25, strand)
            assert not feature.is_contained_by(other)
            assert not other.is_contained_by(feature)
            # b is contained
            other = helpers.DummyFeature(6, 9, strand)
            assert not feature.is_contained_by(other)
            assert other.is_contained_by(feature)
            # only partial overlap
            other = helpers.DummyFeature(6, 19, strand)
            assert not feature.is_contained_by(other)
            assert not other.is_contained_by(feature)
            other = helpers.DummyFeature(1, 7, strand)
            assert not feature.is_contained_by(other)
            assert not other.is_contained_by(feature)
            # edge cases
            other = helpers.DummyFeature(5, 7, strand)
            assert not feature.is_contained_by(other)
            assert other.is_contained_by(feature)
            other = helpers.DummyFeature(7, 10, strand)
            assert not feature.is_contained_by(other)
            assert other.is_contained_by(feature)

    def test_biopython_conversion(self):
        bio = SeqFeature(FeatureLocation(1, 5))
        bio.qualifiers["foo"] = ["bar"]
        # check that features without types are caught
        with self.assertRaises(AssertionError):
            sec = Feature.from_biopython(bio)
        bio.type = "test"
        sec = Feature.from_biopython(bio)
        assert sec.get_qualifier("foo") == tuple(["bar"])
        assert sec.get_qualifier("bar") is None

    def test_string_conversion(self):
        for feature_type in ["cluster", "cds_motif", "test"]:
            for start, end in [(1, 5), (3, 8), (10, 15)]:
                feature = Feature(FeatureLocation(start, end, strand=1),
                                  feature_type=feature_type)
                assert str(feature) == "%s([%d:%d](+))" % (feature_type, start, end)
                feature = Feature(FeatureLocation(start, end, strand=-1),
                                  feature_type=feature_type)
                assert str(feature) == "%s([%d:%d](-))" % (feature_type, start, end)

# since we're about to test assigning to non-slots, shut pylint up
# pylint: disable=assigning-non-slot
    def test_membership(self):
        location = FeatureLocation(0, 1)
        # Features don't have locus tags
        with self.assertRaises(AttributeError):
            Feature(location, feature_type="none").locus_tag = "something"
        # CDSFeatures don't have an 'other_value'
        with self.assertRaises(AttributeError):
            CDSFeature(location, translation="none", gene="a").other_value = 1
        cluster = Cluster(location, 0, 0, products=["a", "b"])
        assert cluster.products == ("a", "b")
        # Clusters have products, not product
        with self.assertRaises(AttributeError):
            cluster.product = ["c", "d"]
# pylint: enable=assigning-non-slot


class TestCluster(unittest.TestCase):
    def create_cluster(self, start, end):
        return Cluster(FeatureLocation(start, end, strand=1),
                       cutoff=1, extent=1, products=['a'])

    def create_cds(self, start, end, strand=1):
        return CDSFeature(FeatureLocation(start, end, strand),
                          locus_tag="%d-%d" % (start, end))

    def setUp(self):
        self.record = Record(Seq("A" * 1000))
        self.start = 100
        self.end = 900
        self.cluster = self.create_cluster(self.start, self.end)
        self.record.add_cluster(self.cluster)
        assert self.cluster.location.start == self.start
        assert self.cluster.location.end == self.end

    def test_trim_unattached(self):
        cluster = self.create_cluster(1, 2)
        cluster.trim_overlapping()
        assert self.cluster.location.start == self.start
        assert self.cluster.location.end == self.end

    def test_trim_empty(self):
        self.cluster.trim_overlapping()
        assert self.cluster.location.start == self.start
        assert self.cluster.location.end == self.end

    def test_trim_contained(self):
        starts = [200, 300, 500]
        ends = [250, 350, 600]
        for start, end in zip(starts, ends):
            feature = self.create_cds(start, end)
            self.record.add_cds_feature(feature)
            assert feature.cluster == self.cluster
        self.cluster.trim_overlapping()
        assert self.cluster.location.start == self.start
        assert self.cluster.location.end == self.end

        for cds in self.record.get_cds_features():
            self.record.remove_cds_feature(cds)
        assert not self.cluster.cds_children

        for end, start in zip(starts, ends):
            feature = self.create_cds(start, end, strand=-1)
            self.record.add_cds_feature(feature)
            assert feature.cluster == self.cluster
        self.cluster.trim_overlapping()
        assert self.cluster.location.start == self.start
        assert self.cluster.location.end == self.end

    def test_trim_leading_overlap(self):
        self.record.add_cds_feature(self.create_cds(self.start - 3, self.start + 3))
        self.record.add_cds_feature(self.create_cds(self.start + 20, self.end - 20))
        self.cluster.trim_overlapping()
        assert self.cluster.location.start == self.start + 3
        assert self.cluster.location.end == self.end

    def test_trim_leading_overlap_with_overlapping_contained(self):  # pylint: disable=invalid-name
        self.record.add_cds_feature(self.create_cds(self.start - 3, self.start + 3))
        self.record.add_cds_feature(self.create_cds(self.start + 1, self.start + 10))
        self.cluster.trim_overlapping()
        assert self.cluster.location.start == self.start + 1
        assert self.cluster.location.end == self.end

    def test_trim_trailing_overlap(self):
        self.record.add_cds_feature(self.create_cds(self.end - 3, self.end + 3))
        self.record.add_cds_feature(self.create_cds(self.start + 20, self.end - 20))
        self.cluster.trim_overlapping()
        assert self.cluster.location.start == self.start
        assert self.cluster.location.end == self.end - 3

    def test_trim_trailing_overlap_with_overlapping_contained(self):  # pylint: disable=invalid-name
        self.record.add_cds_feature(self.create_cds(self.end - 3, self.end + 3))
        self.record.add_cds_feature(self.create_cds(self.end - 10, self.end - 1))
        self.cluster.trim_overlapping()
        assert self.cluster.location.start == self.start
        assert self.cluster.location.end == self.end - 1

    def test_products(self):
        assert self.cluster.products == ("a",)
        self.cluster.add_product("b")
        assert self.cluster.products == ("a", "b")
        with self.assertRaises(AttributeError):
            self.cluster.products.append("c")  # pylint: disable=no-member
        with self.assertRaises(AssertionError):
            self.cluster.add_product(None)
        with self.assertRaises(AssertionError):
            self.cluster.add_product(["C"])


class TestProteinPositionConversion(unittest.TestCase):
    def setUp(self):
        self.func = secmet.feature._convert_protein_position_to_dna

    def test_position_conversion_simple_forward(self):
        location = FeatureLocation(0, 15, strand=1)
        assert len(location) == 15
        assert self.func(0, 2, location) == (0, 6)
        assert self.func(1, 4, location) == (3, 12)

    def test_position_conversion_simple_reverse(self):
        location = FeatureLocation(0, 15, strand=-1)
        assert len(location) == 15
        assert self.func(0, 2, location) == (9, 15)
        assert self.func(1, 4, location) == (3, 12)

    def test_position_conversion_nonzero_start(self):
        location = FeatureLocation(6, 21, strand=1)
        assert len(location) == 15
        assert self.func(0, 2, location) == (6, 12)
        assert self.func(1, 4, location) == (9, 18)

        location = FeatureLocation(6, 21, strand=-1)
        assert len(location) == 15
        assert self.func(0, 2, location) == (15, 21)
        assert self.func(1, 4, location) == (9, 18)

    def test_position_conversion_nonzero_start_compound(self):
        location = CompoundLocation([FeatureLocation(6, 18, strand=1),
                                     FeatureLocation(24, 27, strand=1)])
        assert len(location) == 15
        assert self.func(0, 2, location) == (6, 12)
        assert self.func(1, 4, location) == (9, 18)
        assert self.func(3, 5, location) == (15, 27)

        location = CompoundLocation([FeatureLocation(6, 15, strand=-1),
                                     FeatureLocation(21, 27, strand=-1)])
        assert len(location) == 15
        assert self.func(0, 2, location) == (21, 27)
        assert self.func(1, 4, location) == (9, 24)
        assert self.func(3, 5, location) == (6, 12)


    def test_position_conversion_compound_forward(self):
        location = CompoundLocation([FeatureLocation(0, 6, strand=1),
                                     FeatureLocation(9, 18, strand=1)])
        assert len(location) == 15
        assert self.func(0, 4, location) == (0, 15)
        assert self.func(1, 5, location) == (3, 18)

        location = CompoundLocation([FeatureLocation(0, 6, strand=1),
                                     FeatureLocation(12, 15, strand=1),
                                     FeatureLocation(21, 27, strand=1)])
        assert len(location) == 15
        assert self.func(0, 4, location) == (0, 24)
        assert self.func(1, 5, location) == (3, 27)
        assert self.func(2, 3, location) == (12, 15)

    def test_position_conversion_compound_reverse(self):
        location = CompoundLocation([FeatureLocation(0, 6, strand=-1),
                                     FeatureLocation(9, 18, strand=-1)])
        assert len(location) == 15
        assert self.func(0, 4, location) == (3, 18)
        assert self.func(1, 5, location) == (0, 15)

        location = CompoundLocation([FeatureLocation(0, 6, strand=-1),
                                     FeatureLocation(12, 15, strand=-1),
                                     FeatureLocation(21, 27, strand=-1)])
        assert len(location) == 15
        assert self.func(0, 4, location) == (3, 27)
        assert self.func(1, 5, location) == (0, 24)
        assert self.func(2, 3, location) == (12, 15)

    def test_other(self):
        location = CompoundLocation([FeatureLocation(5922, 6190, strand=1),
                                     FeatureLocation(5741, 5877, strand=1),
                                     FeatureLocation(4952, 5682, strand=1)])
        assert self.func(97, 336, location) == (5243, 6064)

        location = CompoundLocation([FeatureLocation(5922, 6190, strand=-1),
                                     FeatureLocation(5741, 5877, strand=-1),
                                     FeatureLocation(4952, 5682, strand=-1)])
        assert self.func(97, 336, location) == (5078, 5854)

class TestCDSProteinLocation(unittest.TestCase):
    def setUp(self):
        self.magic_split = Seq("ATGGCAxxxxxxGGTxxxxxxATTTGT")
        self.magic = Seq("ATGGCAGGTATTTGT")
        self.translation = "MAGIC"
        self.sub_locations = [FeatureLocation(0, 6, strand=1),
                              FeatureLocation(12, 15, strand=1),
                              FeatureLocation(21, 27, strand=1)]
        self.location = CompoundLocation(self.sub_locations)
        self.cds = CDSFeature(self.location, locus_tag="compound")

    def reverse_strand(self):
        self.magic = self.magic.reverse_complement()
        self.magic_split = self.magic_split.reverse_complement()
        self.sub_locations = [FeatureLocation(loc.start, loc.end, strand=-1) for loc in self.sub_locations]
        self.location = CompoundLocation(self.sub_locations)
        self.cds = CDSFeature(self.location, locus_tag="compound")

    def test_simple_location_forward_complete(self):
        cds = CDSFeature(FeatureLocation(0, 15), locus_tag="simple")
        new = cds.get_sub_location_from_protein_coordinates(0, 5)
        extracted = new.extract(self.magic)
        assert extracted == self.magic
        assert extracted.translate() == self.translation

    def test_simple_location_forward_partial(self):
        cds = CDSFeature(FeatureLocation(0, 15), locus_tag="simple")
        for start, end in [(1, 5), (0, 3), (2, 3), (1, 4)]:
            print("testing", start, end)
            new = cds.get_sub_location_from_protein_coordinates(start, end)
            print(new)
            extracted = new.extract(self.magic)
            assert extracted == self.magic[start * 3: end * 3]
            assert extracted.translate() == self.translation[start:end]

    def test_compound_location_forward_full(self):
        new = self.cds.get_sub_location_from_protein_coordinates(0, 5)
        assert isinstance(new, CompoundLocation)
        assert len(new.parts) == 3
        print(list(map(str, self.cds.location.parts)))
        print(list(map(str, new.parts)))
        assert len(new) == len(self.cds.location)
        assert new == self.location, "%s != %s" % (str(new), str(self.location))
        extracted = new.extract(self.magic_split)
        assert extracted == self.magic
        assert extracted.translate() == self.translation[0:5]

    def test_compound_forward_within_single(self):
        new = self.cds.get_sub_location_from_protein_coordinates(0, 2)
        assert isinstance(new, FeatureLocation)
        assert len(new) == 6
        assert new.start == 0
        assert new.end == 6
        assert new.extract(self.magic_split).translate() == self.translation[0:2]

        new = self.cds.get_sub_location_from_protein_coordinates(2, 3)
        assert isinstance(new, FeatureLocation)
        assert len(new) == 3
        assert new.start == 12
        assert new.end == 15
        assert new.extract(self.magic_split).translate() == self.translation[2:3]

    def test_compound_forward_over_multiple(self):
        new = self.cds.get_sub_location_from_protein_coordinates(2, 4)
        assert isinstance(new, CompoundLocation)
        print(list(map(str, self.cds.location.parts)))
        print(list(map(str, new.parts)))
        assert len(new.parts) == 2
        assert len(new) == 6
        assert new.parts[0].start == 12
        assert new.parts[0].end == 15
        assert new.parts[1].start == 21
        assert new.parts[1].end == 24
        assert new.extract(self.magic_split).translate() == self.translation[2:4]

    def test_compound_location_reverse_full(self):
        self.reverse_strand()
        for direction in [1, -1]:
            print("sub locations sorted in direction:", direction)
            cds = CDSFeature(CompoundLocation(self.sub_locations[::direction]), locus_tag="compound")
            new = cds.get_sub_location_from_protein_coordinates(0, 5)
            assert isinstance(new, CompoundLocation)
            assert len(new.parts) == 3
            print(list(map(str, cds.location.parts)))
            print(list(map(str, new.parts)))
            assert len(new) == len(cds.location)
            assert new.extract(self.magic_split).translate() == self.translation[0:5]

    def test_compound_location_reverse_single(self):
        self.reverse_strand()
        for direction in [1, -1]:
            print("sub locations sorted in direction:", direction)
            cds = CDSFeature(CompoundLocation(self.sub_locations[::direction]), locus_tag="compound")

            new = cds.get_sub_location_from_protein_coordinates(0, 2)
            assert isinstance(new, FeatureLocation)
            assert len(new) == 6
            assert new.start == 21
            assert new.end == 27
            assert new.extract(self.magic_split).translate() == self.translation[0:2]

            new = cds.get_sub_location_from_protein_coordinates(2, 3)
            assert isinstance(new, FeatureLocation)
            assert len(new) == 3
            assert new.start == 12
            assert new.end == 15
            assert new.extract(self.magic_split).translate() == self.translation[2:3]

    def test_compound_location_reverse_multiple(self):
        self.reverse_strand()
        for direction in [1, -1]:
            print("sub locations sorted in direction:", direction)
            cds = CDSFeature(CompoundLocation(self.sub_locations[::direction]), locus_tag="compound")

            new = cds.get_sub_location_from_protein_coordinates(2, 4)
            assert isinstance(new, CompoundLocation)
            print(list(map(str, cds.location.parts)))
            print(list(map(str, new.parts)))
            assert len(new.parts) == 2
            assert len(new) == 6
            assert new.parts[0].start == 12
            assert new.parts[0].end == 15
            assert new.parts[1].start == 3
            assert new.parts[1].end == 6
            assert new.extract(self.magic_split).translate() == self.translation[2:4]

    def test_complicated(self):
        parts = [FeatureLocation(121124, 122061, 1), FeatureLocation(122339, 122383, 1),
                 FeatureLocation(122559, 122666, 1), FeatureLocation(122712, 122874, 1),
                 FeatureLocation(123060, 123337, 1), FeatureLocation(123481, 123749, 1),
                 FeatureLocation(123809, 124032, 1), FeatureLocation(124091, 124193, 1),
                 FeatureLocation(124236, 124401, 1), FeatureLocation(124684, 124724, 1)]
        location = CompoundLocation(parts, operator="join")
        cds = CDSFeature(location, locus_tag="complicated")
        seq = ("ATGAGCCCTCGTCTAGACTACAATGAAGGATACGATTCCGAAGACGAGGAGATCCCCCGTTACGTACACCAT"
               "TCTAGAGGAAAGAGTCATAGATCCGTGAGGACGTCAGGTCGCTCACGCACGTTGGATTACGACGGGGATGAT"
               "GAAGCTAGTGACCACGCTGCCCCCTCCGGGATTGATCGGGACGCTCGAGCCTGTCCAACATCTCGCAGATAT"
               "ACTGATGACTGCCTTGAGACACATAAATTTCGAGGTGCCCGCTCCTCTCGCTCCCGTGGACGAACCGATGAT"
               "AACAAGGTTTTGTACTACACCAAGTATCGCAGCCCGGCTAAGGACTTGCCTATCGAGCGTGATCCCGAGGGT"
               "ATTAATTTATTCAAGGTCCGACAGCACACACGGCCAAGTGACGCTCATGTGCCCAGTGGATACCGTGAGCCC"
               "TACGAAGTCAAGGTCGACGAGTATGAGGATGATCATCCCCGTACATGCACTAGCCGCCGTGACTCTAGACAG"
               "CCGAAAGTCTACAAGGTCCGGGTTGATGAGTACGAGGATAACCTCCCTGCACGCTCTCACACTGACTTTCGC"
               "GAGTCTCCACGGTCTGAAAGATGCTCTAGCCGCTACACCGAGGACTCGAAGCCTGGGGAGCTTCCTCCCCGC"
               "TCAGGGCCCTGTCGGTCCAGCAGGCCTTCTCCGGTCGATGAGGACGTCGAGTATGAGATCCGTGAGCCCCGA"
               "GGGCATCGCTCCAGTCGACACTCTACAGATGTTGACTTTCAGCCAGTAGAACAACATCCTCGCTTTGGACAA"
               "CGTGGACTCAGCAGACCTTCGCGGGTTGATGAGGAAGTCGATTATGAGATCCGTGAGCCCCGTGGCAATCGT"
               "GTCAGTCACGCTGCTCATGGTGACAGCCCCTGTCAGGACCAAAGCTCCAGGCATATCGGCATTCAATTGTGG"
               "AGTACGCGCGGACCCCGGGCGGCTGGCCGTGGCCGGGGTCCTGATGAGTCTGACGATGTTGAGCCCTAGGCA"
               "GGGAATTGCCGTAATGCTCTTCAAACTGTATAGCAAGCTCAGCATCAATTCTTTAACTGGCAGGCGCTCTGC"
               "TCGCGCGTTTCTCTCTTGGGGTGGTTGGTTTGACTGTAGATTTCCTCTTTCAAGGCTTCTAGATACACCTTT"
               "GGAAGATAGCAACGCTATGCAAGATATTTTTGATAATTCAAATCCTTTTTACACATGGAATAGCTGGTGTTC"
               "CTGTTTTATCTAGGCAATTGACCCACGCCATCTCGGTAGGTACGGTAAAAGCAAGCCGTAATCTCGTATGGC"
               "TTCATCCTTAGCATCGTATAGATCTCCACTCGGGACTCGGCCAGGGATCTTCCATCAATCAACGTGAAGAAG"
               "TCCAGCACCCCGCTGAATCATAATATCCTACCGATTCTGCTCTCTTCACCTCTAGATACCCCTCTAGACTCC"
               "TGTCAACATGTTCCGTACAGTCGAAGACCGCCCGACCCCAAAAGAGGTATATAACTGGCGGCTGTACACCGA"
               "GGCCACCATCATTGCCACTGGTACACTCTTGTGAGTAGGTGCTGTTGTAACGAAAAACATCCAACTGATCCG"
               "CCAGGTTCGGCTATGACTCGGCTTTTGTGGGAACTACCATTGCCCGCCAAAGCTTCGTTGATGCCTTCAACA"
               "TCGTCGAGTCGGAGGCGGCGGATATTTCAAGCAATATCACGTCAACCTTTCAGGCCGGCGCATTTTTCGGCG"
               "CCATCTTCTGCTTCTTGCCTGAGTGAAGCCGTTAGAGACGGTCTCACTGGCTAACCGGACCAAGTGACCGAC"
               "AAAATTGGGCGTAAATGGGCCCTTCAGGCAAACACACTGCTGTTTCTTATTGGCGCGATTGTGATGACGGCT"
               "GCAACACATCACCTTTCCTATATATGTAAGTCATATCCCCGTAGTAGTCAAGGTTGTTAACTAGAGCAGATG"
               "CTGGACGAGCTCTCACCGGCATCGCATGCGGCGCTATCACCGCGACCGTCCCCAGCTATATTGCCGAGCTGT"
               "CAATCGTGTCGATCCGGGGCTTCCTCACCGGGTTCTTCGAAGTCGCATACCAGATTGGTAGCTTGGTTGGAT"
               "TCTGGATCAACTATGGCATTAACGAGAACATGGACAACTCCTCGGCCGCAAGCTGGAGAGTGCCTATGGCAG"
               "TCCAGATCATCCCCGCAGGAGTCCTTTTCATTGGTGGCTTTTCCTCCATGAGAGTCCTCTCTGGCTGATGCG"
               "AAAAGACAGTGAGGATGCCGCGACGGCTGCCCTGGAGGCGTTGAGGAAACTGCCACGGTCTCATCAATGTAA"
               "TCTCCCACCAAGACTCAGGACATAGTCCCATGCTGACTATTTTAGATGTCCAGGAAGACATCGAGATGAACC"
               "GCACCAGGCTGCTGGAGGAAGCTCGGATCGCCGAGAAGTACGGACAAGGTTGGTTGGCATATATCCGAGGCG"
               "CACTCTTCGAGCTCTCGCGCCATGGGATGTGGAATCGTGTTCTGCTCGTCCTCTGTGCCTTTGCACTGCAGA"
               "ATATGTCGGGAGCTGCTGCTATCAACTACTATTCCCCCATACTCTTTGCGTCGTTGGGGATCACTGATGTCG"
               "CTCTGTATACAGGTATTTATGGCCTGGTAAAAGGTAAGTTCTTCTCCTTAAGTATCTCTGGCTGACAATAGG"
               "GATTAACTGATGAGTTTACAGCCGTCGCATCAATTATATTCTACGGCATTCTCATTGATATGTGGGGCCGCC"
               "GACGTCCGACCATTGTTTCGTCACTGGCCTGCCCTCTATGTCTCTGGTTTGTGGGTGCATACGTCAAAGTTG"
               "GGCATCCAGCCGATATCATAGACGCCGGCGGGGAATTGTCCCCCTCCACGGAGGCTGGTGGTAGAGCGGCGA"
               "CTGCGATGATTATGATCTACTCCGTCTTGTAAGTGCCCCTCACTTTTGAATGGGCTTCAGCTTGGAACTCGA"
               "GTAACTGGTATCCAGTTGGTCTTTTGGTCTCAACGGTATCCCCTGGATTGTCTCCGCCGAAATCTTCCCCGG"
               "CGCGCTGCGAAATCTCACGGGGACATGGGCTGCGCTGGTGCAATGGTATGCAATTCCCTTCACCTAGTATCC"
               "ATATCTAAATCAGCAGGTTGATCCAATTCGTTATCACCAAAGCTCTCCCGTACATCTTCAATAGCCTTGGGT"
               "ACGGGACGTGGTTCTTCTTCGCCTCCTGGATGCTGCTCGCTATCATTTGGTCATTCTTTTTTCTCCCGGAAA"
               "CCAAGGGGAAGACTCTCGATGAAATGCATACGATCTTGTACGTTTCTCTCCGTCGAAATGTGGTCTTGGCTA"
               "ATGAATCAGCGGCCATTCTCTCGCCGAAGAGCAGGGTAAGGGTGAGGTTCGAGATAACACTACTAAAAGTGA"
               "TCGGGAGGCTGTCTAGTCCAGTAGTTCTAGAGGACTATTGGCTGGATGATTCCTCTGATGATTTTTGATTGG"
               "TGGTGAAAATGTTGGATGTTTAATGCCAATGTACTGGGAGAGAACATGCCGATAGTACATACCGCTGTGTTG"
               "TATATCGAAGACGGTTGATTTATATATCTTAGTCTTTCAAAAGACGGCACTCACACAATCACACTTCGATGA")
        translation = ("MSPRLDYNEGYDSEDEEIPRYVHHSRGKSHRSVRTSGRSRTLDYDGDDEASDHAAPSGIDRDAR"
                       "ACPTSRRYTDDCLETHKFRGARSSRSRGRTDDNKVLYYTKYRSPAKDLPIERDPEGINLFKVRQ"
                       "HTRPSDAHVPSGYREPYEVKVDEYEDDHPRTCTSRRDSRQPKVYKVRVDEYEDNLPARSHTDFR"
                       "ESPRSERCSSRYTEDSKPGELPPRSGPCRSSRPSPVDEDVEYEIREPRGHRSSRHSTDVDFQPV"
                       "EQHPRFGQRGLSRPSRVDEEVDYEIREPRGNRVSHAAHGDSPCQDQSSRHIGIQLWTGVPVLSR"
                       "QLTHAISTPVNMFRTVEDRPTPKEVYNWRLYTEATIIATGTLLFGYDSAFVGTTIARQSFVDAF"
                       "NIVESEAADISSNITSTFQAGAFFGAIFCFLPEADAGRALTGIACGAITATVPSYIAELSIVSI"
                       "RGFLTGFFEVAYQIGSLVGFWINYGINENMDNSSAASWRVPMAVQIIPAGVLFIGGFSSMREDI"
                       "EMNRTRLLEEARIAEKYGQGWLAYIRGALFELSRHGMWNRVLLVLCAFALQNMSGAAAINYYSP"
                       "ILFASLGITDVALYTGIYGLVKAVASIIFYGILIDMWGRRRPTIVSSLACPLCLWFVGAYVKVG"
                       "HPADIIDAGGELSPSTEAGGRAATAMIMIYSVFWSFGLNGIPWIVSAEIFPGALRNLTGTWAAL"
                       "VQWLIQFVITKALPYIFNSLGYGTWFFFASWMLLAIIWSFFFLPETKGKTLDEMHTIFLSKDGT"
                       "HTITLR")
        new = cds.get_sub_location_from_protein_coordinates(353, 412)
        # pad the beginning to match the location
        assert new.extract(Seq("x" * location.start + seq)).translate() == translation[353:412]


class TestCDSFeature(unittest.TestCase):
    def test_required_identifiers(self):
        with self.assertRaises(ValueError):
            CDSFeature(FeatureLocation(1, 5))
        assert CDSFeature(FeatureLocation(1, 5), locus_tag="foo")
        assert CDSFeature(FeatureLocation(1, 5), protein_id="foo")
        assert CDSFeature(FeatureLocation(1, 5), gene="foo")


class TestGeneFunction(unittest.TestCase):
    def test_membership(self):
        assert GeneFunction.OTHER
        with self.assertRaises(AttributeError):
            print(GeneFunction.non_existant)

    def test_equality(self):
        assert GeneFunction.OTHER == GeneFunction.OTHER
        assert GeneFunction.CORE != GeneFunction.OTHER

    def test_string_conversion(self):
        assert str(GeneFunction.CORE) == "biosynthetic"
        assert str(GeneFunction.ADDITIONAL) == "biosynthetic-additional"
        assert str(GeneFunction.OTHER) == "other"
        for member in dir(GeneFunction):
            if not member.isupper() or member in ["CORE", "ADDITIONAL"]:
                continue
            assert str(getattr(GeneFunction, member)) == member.lower()

    def test_cds_function(self):
        cds = CDSFeature(FeatureLocation(1, 5), locus_tag="foo")
        # default value
        assert cds.gene_functions.get_classification() == GeneFunction.OTHER
        assert cds.gene_function == GeneFunction.OTHER
        # check bad values can't be assigned
        with self.assertRaises(AssertionError):
            cds.gene_functions.add("other", "a", "b")
        with self.assertRaises(AttributeError):
            cds.gene_functions = 0

        cds.gene_functions.add(GeneFunction.ADDITIONAL, "first_tool", "dummy")
        assert cds.gene_functions.get_classification() == GeneFunction.ADDITIONAL
        assert cds.gene_function == GeneFunction.ADDITIONAL
        # conflicting, so back to OTHER
        cds.gene_functions.add(GeneFunction.TRANSPORT, "other_tool", "dummy")
        assert cds.gene_functions.get_classification() == GeneFunction.OTHER
        assert cds.gene_function == GeneFunction.OTHER
        # but smcogs overrides that
        cds.gene_functions.add(GeneFunction.REGULATORY, "smcogs", "dummy")
        assert cds.gene_functions.get_classification() == GeneFunction.REGULATORY
        # and cluster definition overrides even that
        cds.gene_functions.add(GeneFunction.CORE, "cluster_definition", "dummy")
        assert cds.gene_functions.get_classification() == GeneFunction.CORE

        # and that we still have tracked these
        smcogs = cds.gene_functions.get_by_tool("smcogs")
        assert len(smcogs) == 1
        assert smcogs[0].function == GeneFunction.REGULATORY

        adds = cds.gene_functions.get_by_function(GeneFunction.ADDITIONAL)
        assert len(adds) == 1
        assert adds[0].tool == "first_tool"

    def test_cds_function_conversion(self):
        cds = CDSFeature(FeatureLocation(1, 5), locus_tag="foo")
        assert cds.gene_function == GeneFunction.OTHER
        assert CDSFeature.from_biopython(cds.to_biopython()[0]).gene_function == GeneFunction.OTHER
        cds.gene_functions.add(GeneFunction.ADDITIONAL, "tool", "desc")
        assert cds.gene_function == GeneFunction.ADDITIONAL
        assert CDSFeature.from_biopython(cds.to_biopython()[0]).gene_function == GeneFunction.ADDITIONAL
