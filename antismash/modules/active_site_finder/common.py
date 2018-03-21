# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import re
from typing import List, Optional

from antismash.common import fasta, path, secmet, subprocessing, utils


class Alignment:
    def __init__(self, domain: secmet.feature.Domain, query: str, hit: str, hit_start: int, hit_end: int):
        assert isinstance(domain, secmet.feature.Domain)
        self.domain = domain
        self.query = str(query)
        self.hit = str(hit)
        self.hit_start = int(hit_start)
        self.hit_end = int(hit_end)

    def get_signature(self, positions: List[int], expected: List[str]) -> str:
        assert len(positions) == len(expected)
        return get_signature(self.query, self.hit, positions, expected)

    def extract_position(self, position: int) -> str:
        assert 0 <= self.hit_start < position < self.hit_end
        return self.hit[position - self.hit_start - 1]

    def extract_positions(self, positions: List[int]) -> str:
        return "".join(self.extract_position(position) for position in positions)


class ActiveSiteAnalysis:
    def __init__(self, target_domain: str, candidates: List[secmet.feature.Domain],
                 database: str, command: str, positions: List[int],
                 expected_values: List[str], emissions: Optional[List[float]] = None) -> None:
        self.target_domain = str(target_domain)
        self.database = path.get_full_path(__file__, 'data', database)
        self.command = str(command)
        assert self.command in ["hmmpfam2", "hmmscan"]
        self.positions = list(map(int, positions))
        self.expected_values = list(map(str, expected_values))
        if len(expected_values) != len(positions):
            raise ValueError("Number of expected values must match number of positions")
        self.emissions = None
        if emissions:
            self.emissions = list(map(int, emissions))
            if len(self.emissions) != len(positions):
                raise ValueError("Number of emissions must match number of positions")

        self.domains_of_interest = []
        for candidate in candidates:
            if not isinstance(candidate, secmet.feature.Domain):
                raise TypeError("Candidates must be domains, not %s" % type(candidate))
            if candidate.domain == self.target_domain:
                self.domains_of_interest.append(candidate)

    def get_alignments(self) -> List[Alignment]:
        if not self.domains_of_interest:
            return []


        # TODO
        # run_hmmpfam fails for some cases due to a biopython bug.
        # run each domain on its own and warn about the bad ones, then continue



        # for safety of the tools, rename long domain names to a simple numeric index
        data = fasta.get_fasta_from_features(self.domains_of_interest, numeric_names=True)
        assert data, "empty fasta created"
        if self.command == "hmmpfam2":
            extra_args = ["-T", "0",  # min score
                          "-E", "0.1"]  # max evalue
            logging.critical("running hmmpfam2 -T 0 -E 0.1 %s -\n %s", self.database, data)
            try:
                results = subprocessing.run_hmmpfam2(self.database, data, extra_args=extra_args)
            except:
                print(self.domains_of_interest[0])
                raise
        else:  # hmmscan
            extra_args = ["-domE", "0.1"]  # max evalue
            results = subprocessing.run_hmmscan(self.database, data, opts=extra_args)

        alignments = []
        for result in results:
            if not result.hsps:
                continue
            assert result.id == result.hsps[0].aln[0].id
            if not result.hsps:
                continue
            # fetch back the real domain from the numeric index used in the fasta
            domain = self.domains_of_interest[int(result.id)]
            alignments.append(Alignment(domain, result.hsps[0].aln[0].seq, result.hsps[0].aln[1].seq,
                                        result.hsps[0].hit_start, result.hsps[0].hit_end))
        return alignments

    def scaffold_matches(self, alignment: Alignment):
        return get_signature(alignment.query, alignment.hit, self.positions, self.expected_values) == "".join(self.expected_values)


def get_signature(query, hmm, positions, expected=None) -> str:
    """

        positions are 1-indexed
    """
    ungapped = str(hmm).replace('.', '')
    if positions[-1] > len(ungapped):
        logging.warning("scaffold coordinate %s outside hsp!", positions[-1] - 1)
        return None
    ref_signature = "".join(ungapped[pos - 1] for pos in positions)
    if expected:
        assert ref_signature == "".join(expected)
    signature = utils.extract_by_reference_positions(query, hmm, [pos - 1 for pos in positions])
    return signature


def get_scaffold_annotation(result, positions, values, expected_emissions=None):
    """ Generate annotation from scaffold information

        positions are 1-indexed
    """

    query_seq = result.hsps[0].aln[0].seq
    hmm_seq = result.hsps[0].aln[1].seq
    assert len(values) == len(positions)
    if not expected_emissions:
        emissions = ["n.d."] * len(positions)
    else:
        emissions = list(map(int, expected_emissions))  # TODO: should floats (0., 1.) really be converted to ints?
    assert len(emissions) == len(positions)

    # Check scaffold matches
    matches = []

    signature = get_signature(query_seq, hmm_seq, positions, values)
    assert len(signature) == len(positions), "%d != %d" % (len(signature), len(positions))
    for i, char in enumerate(signature):
        position = positions[i] - 1
        expected = values[i]

        # We have to use a RegEx here to allow negations and more complex queries; ignore case (?i)
        match = bool(re.match("(?i)" + expected, char))
        matches.append(str(match))
        logging.debug("Scaffold coordinate %s; query aa %s; "
                      "scaffold value %s; emission probability %s; match %s",
                      position, char, expected, emissions[i], match)

    overall_match = all(matches)

    logging.debug("Overall Scaffold Match: %s", str(overall_match).upper())

    # Generate Feature qualifiers
    return ("Scaffold coordinates: (%s); scaffold residues: (%s); expected: (%s); matchArray: (%s); "
            "emission probability array (%s); overall match: %s") % (
            ",".join(map(str, positions)), signature, ",".join(values),
            ",".join(map(str, matches)), ",".join(map(str, emissions)),
            str(overall_match).upper())


def get_prediction_annotation(result, positions, values, result_label, comment, expected_emissions=None):
    "gererate annotation from choices/prediction information, a single 'choice' at a time"

    query_seq = result.hsps[0].aln[0].seq
    hmm_seq = result.hsps[0].aln[1].seq
    if not expected_emissions:
        emissions = ["n.d."] * len(positions)
    else:
        emissions = list(map(int, expected_emissions))  # TODO: again, should floats (0., 1.) really be converted to ints?

    signature = get_signature(query_seq, hmm_seq, positions, values)

    logging.debug("testing %s (%s):", result_label, comment)

    matches = []
    for i, char in enumerate(signature):
        value = values[i]

        position = int(positions[i]) - 1

        matches.append(bool(re.match("(?i)" + value, char)))

        logging.debug("Offset %s; Expected: %s; observed in query %s; Emission %s; match %s",
                      position,
                      value,
                      char,
                      emissions[i],
                      matches[-1])

    overall_match = all(matches)

    logging.debug("Overall Match for prediction %s: %s", result_label, str(overall_match).upper())
    logging.debug("================================")

    description = (("Description: %s, choice result: %s, choice coordinates: (%s); residues: (%s); "
                   "expected for choice: (%s); matchArray: (%s); emission probability array (%s); overall match: %s") %
                   (comment,
                    result_label,
                    ",".join(map(str, positions)),
                    signature,
                    ",".join(values),
                    ",".join(map(str, matches)),
                    ",".join(map(str, emissions)),
                    str(overall_match).upper()))

    choice_string = ""
    if overall_match:
        choice_string = "Full match for prediction: %s" % result_label
    return description, choice_string
