# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import re
import os
from typing import List

import Bio
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from helperlibs.bio import seqio

from antismash.common import gff_parser
from antismash.common.secmet import Record, CDSFeature
from antismash.config import update_config

from .utils import generate_unique_id

def parse_input_sequence(filename, minimum_length=-1, start=-1, end=-1) -> List[Record]:
    """ Parse input records contained in a file

        Arguments:
            filename: the path of the file to read
            minimum_length: records with length less than this will be ignored
                            if not positive, all records are included

        Returns:
            A list of secmet.Record instances, one for each record in the file
    """
    logging.info('Parsing input sequence %r', filename)
    if not isinstance(minimum_length, int):
        raise TypeError("minimum_length must be an int")

    records = []
    if not os.path.exists(filename):
        msg = "Sequence file not found: %r" % filename
        logging.error(msg)
        raise ValueError(msg)

    try:
        record_list = list(seqio.parse(filename))
        if not record_list:
            raise RuntimeError('No records could be read from file %r' % filename)
        for record in record_list:
            if minimum_length < 1 \
                    or len(record.seq) >= minimum_length \
                    or 'contig' in record.annotations \
                    or 'wgs_scafld' in record.annotations \
                    or 'wgs' in record.annotations:
                records.append(record)
    except (ValueError, AssertionError) as err:
        logging.error('Parsing %r failed: %s', filename, err)
        raise
    except Exception as err:
        logging.error('Parsing %r failed with unhandled exception: %s',
                      filename, err)
        raise

    # before conversion to secmet records, trim if required
    if start > -1 or end > -1:
        if len(records) > 1:
            raise ValueError("--start and --end options cannot be used with multiple records")
        trim_sequence(records[0], start, end)

    return [Record.from_biopython(record) for record in records]

def pre_process_sequences(sequences, options, genefinding) -> List[Record]:
    """ hmm

        - gaps removed
        - record ids adjusted to be unique
        - record ids are valid

        Note: Record instances will be altered in-place.

        Arguments:
            sequences: the secmet.Record instances to process
            options: an antismash Config instance
            genefinding: the module to use for genefinding, must have
                         run_on_record() implemented

        Returns:
            A list of altered secmet.Record
    """
    # keep count of how many records matched filter
    matching_filter = 0

    # Check if seq_records have appropriate content
    for i, sequence in enumerate(sequences):
        if options.limit_to_record and options.limit_to_record != sequence.id:
            sequence.skip = "did not match filter: %s" % options.limit_to_record
        else:
            matching_filter += 1

        sequence.record_index = i
        sequence.seq = Seq(str(sequence.seq).replace("-", "").replace(":", ""))
        # Check if seq_record has either a sequence or has at least 80% of CDS features with 'translation' qualifier
        cdsfeatures = sequence.get_cds_features()
        cdsfeatures_with_translations = sum([1 for cdsfeature in cdsfeatures if cdsfeature.translation])
        assert cdsfeatures_with_translations == len(cdsfeatures)
        if not sequence.seq or (options.input_type == 'nucl' and \
                                not str(sequence.seq).replace("N", "")):
            logging.error("Record %s has no sequence, skipping.", sequence.id)
            sequence.skip = "contains no sequence"
            continue

        if options.input_type == 'prot':
            if is_nucl_seq(sequence.seq):
                logging.error("Record %s is a nucleotide record, skipping.", sequence.id)
                sequence.skip = "nucleotide record in protein mode"
                continue
        elif options.input_type == 'nucl':
            if not isinstance(sequence.seq.alphabet, Bio.Alphabet.NucleotideAlphabet) and not is_nucl_seq(sequence.seq):
                logging.error("Record %s is a protein record, skipping.", sequence.id)
                sequence.skip = "protein record in nucleotide mode"
                continue
            sequence.seq.alphabet = Bio.Alphabet.generic_dna

    if options.limit_to_record:
        limit = options.limit_to_record
        if matching_filter == 0:
            logging.error("No sequences matched filter: %s", limit)
            raise ValueError("No sequences matched filter: %s" % limit)
        elif matching_filter != len(sequences):
            logging.info("Skipped %d sequences not matching filter: %s",
                         len(sequences) - matching_filter, limit)

    #If protein input, convert all protein seq_records to one nucleotide seq_record
    if options.input_type == 'prot':
        sequences = generate_nucl_seq_record(sequences)

    # catch WGS master or supercontig entries
    if records_contain_shotgun_scaffolds(sequences):
        raise RuntimeError("Incomplete whole genome shotgun records are not supported")

    # Now remove small contigs < minimum length again
    for sequence in sequences:
        if len(sequence.seq) < options.minlength:
            sequence.skip = "smaller than minimum length (%d)" % options.minlength

    # Make sure we don't waste weeks of runtime on huge records, unless requested by the user
    warned = False
    if options.limit > -1:
        meaningful = 0
        for sequence in sequences:
            if sequence.skip:
                continue
            meaningful += 1
            if meaningful > options.limit:
                if not warned:
                    logging.warning("Only analysing the first %d records (increase via --limit)", options.limit)
                    warned = True
                sequence.skip = "skipping all but first {0} meaningful records (--limit {0}) ".format(options.limit)

    options = update_config({"triggered_limit" : warned}) # TODO is there a better way

    # Check if no duplicate locus tags / gene IDs are found
    ensure_no_duplicate_gene_ids(sequences)

    # Check GFF suitability
    if options.genefinding_gff3:
        single_entry = gff_parser.check_gff_suitability(options, sequences)

    all_record_ids = {seq.id for seq in sequences}
    # Ensure all records have unique names
    if len(all_record_ids) < len(sequences):
        all_record_ids = set()
        for record in sequences:
            if record.id in all_record_ids:
                record.id = generate_unique_id(record.id, all_record_ids)[0]
            all_record_ids.add(record.id)
        assert len(all_record_ids) == len(sequences), "%d != %d" % (len(all_record_ids), len(sequences))
    # Ensure all records have valid names
    for record in sequences:
        fix_record_name_id(record, all_record_ids)

    for sequence in sequences:
        if sequence.skip:
            continue
        if len(sequence.get_cds_features()) < 1:
            if options.genefinding_gff3:
                logging.info("No CDS features found in record %r but GFF3 file provided, running GFF parser.", sequence.id)
                gff_parser.run(sequence, single_entry, options)
                # since gff_parser adds features, make sure they don't share ids
                ensure_no_duplicate_gene_ids(sequences)
            elif options.genefinding_tool != "none":
                logging.info("No CDS features found in record %r, running gene finding.", sequence.id)
                genefinding.run_on_record(sequence, options)
            if len(sequence.get_cds_features()) < 1:
                logging.info("No genes found, skipping record")
                sequence.skip = "No genes found"
                continue
        #Fix locus tags
        fix_locus_tags(sequence)

    ensure_gap_notation_consistent(sequences)

    return sequences


def ensure_gap_notation_consistent(records) -> None:
    """ Ensures all sequences use N for gaps instead of -, X, or x

        Arguments:
            records: the secmet.Records to alter

        Returns:
            None
    """
    for record in records:
        for char in "-Xx":
            new_seq = str(record.seq).replace(char, "N")
            record.seq = Seq(new_seq, alphabet=record.seq.alphabet)


def trim_sequence(record, start, end) -> SeqRecord:
    """ Trims a record to the range given

        Arguments:
            record: the Bio.SeqRecord to trim
            start: the start position (inclusive)
            end: the end position (exclusive)

        Returns:
            A new, shortened Bio.SeqRecord instance
    """
    if start >= len(record):
        raise ValueError('Specified analysis start point of %r is outside record' % start)
    if end > len(record):
        raise ValueError('Specified analysis end point of %r is outside record' % end)
    if end > -1 and end <= start:
        raise ValueError("Trim region start cannot be greater than or equal to end")

    if start < 0:
        start = 0
    if end < 0:
        end = len(record)
    return record[start:end]


def is_nucl_seq(sequence) -> bool:
    """ Determines if a sequence is a nucleotide sequence based on content.

        Arguments:
            sequence: the sequence to check, either a string or Bio.Seq

        Returns:
            True if less than 20% of bases are not a,c,g,t or n
    """
    other = str(sequence).lower()
    for char in "acgtn":
        other = other.replace(char, "")
    return len(other) < 0.2 * len(sequence)


def generate_nucl_seq_record(records) -> Record:
    """ Generate single nucleotide-based record from multiple records

        Arguments:
            records: the records from which to generate the record

        Returns:
            a new secmet.Record instance containing all record info
    """
    if not records:
        raise ValueError("Cannot generate nucleotide records of empty input")

    record = Record(Seq(""), id="Protein_Input", name="ProteinInput", # TODO: surely we can use a better base id/name/description
                   description="antiSMASH protein input")
    position = 0
    cds_names = set()
    for record in records:
        start = position
        end = position + len(record) * 3
        position += len(record) * 3 + 1000
        location = FeatureLocation(start, end, strand=1)
        name = record.id[:15].replace(" ", "_")
        if name in cds_names:
            name, _ = generate_unique_id(name[:8], cds_names)
        cds_names.add(name)
        translation = str(record.seq).replace('.', 'X')
        cdsfeature = CDSFeature(location, translation, product=name, locus_tag=name)
        record.add_cds_feature(cdsfeature)
    return record


def records_contain_shotgun_scaffolds(records) -> bool:
    """ Check if given records contain a WGS master record or supercontig record

        Arguments:
            records: an iterable of secmet.Record

        Returns:
            True if one of the given records is a WGS or supercontig record
    """
    for record in records:
        if 'wgs_scafld' in record.annotations \
                or 'wgs' in record.annotations \
                or 'contig' in record.annotations:
            return True
    return False


def ensure_no_duplicate_gene_ids(sequences) -> None:
    """ Ensures that every CDS across all sequences has a unique id

        Arguments:
            sequences: the secmet.Record instances to process

        Returns:
            None
    """
    no_tag = "no_tag_found"

    high_water_mark = 0
    all_ids = set()
    for sequence in sequences:
        for cdsfeature in sequence.get_cds_features():
            name = cdsfeature.get_name()
            if not name:
                name = no_tag
            if name == no_tag:
                name, high_water_mark = generate_unique_id(name[:8], all_ids,
                                                start=high_water_mark + 1)
            elif name in all_ids:
                name, _ = generate_unique_id(name[:8], all_ids, start=1)
            if cdsfeature.product is None:
                cdsfeature.product = name
            cdsfeature.locus_tag = name
            all_ids.add(name)


def fix_record_name_id(record, all_record_ids) -> None:
    """ Changes a record's name and id to be no more than 16 characters long,
        so it can be used in GenBank files.

        If record name is too long, the prefix c000X is used

        Arguments:
            record: the record to alter
            all_record_ids: a set of all known record ids

        Returns:
            None
    """

    def _shorten_ids(idstring):
        contigstrmatch = re.search(r"onti?g?(\d+)\b", idstring)
        if not contigstrmatch:
            # if there is a substring "[Ss]caf(fold)XXX" use this number
            contigstrmatch = re.search(r"caff?o?l?d?(\d+)\b", idstring)
        if not contigstrmatch:
            # if there is a substring " cXXX" use this number
            contigstrmatch = re.search(r"\bc(\d+)\b", idstring)
        if contigstrmatch:
            contig_no = int(contigstrmatch.group(1))
        else:
            # if the contig number cannot be parsed out, just count the contigs from 1 to n
            contig_no = record.record_index

        return "c{ctg:05d}_{origid}..".format(ctg=contig_no, origid=idstring[:7])

    old_id = record.id

    if record.id in ["unknown", "unknown.1"]:  #TODO oddly specific
        old_id = record.id
        record.id = "unk_seq_{ctg:05d}".format(ctg=record.record_index)
        logging.warning('Invalid sequence id "%s", replaced by %s', old_id, record.id)

    if len(record.id) > 16:
        # Check if it is a RefSeq accession number like NZ_AMZN01000079.1 that
        # is too long just because of the version number behind the dot
        if (record.id[-2] == "." and
                record.id.count(".") == 1 and
                len(record.id.partition(".")[0]) <= 16 and
                record.id.partition(".")[0] not in all_record_ids):
            record.id = record.id.partition(".")[0]
            all_record_ids.add(record.id)
        else: #Check if the ID suggested by _shorten_ids is unique
            if _shorten_ids(old_id) not in all_record_ids:
                name = _shorten_ids(old_id)
            else:
                name, _ = generate_unique_id(record.id[:12], all_record_ids, max_length=16)
            record.id = name
            all_record_ids.add(name)

        logging.warning('Fasta header too long: renamed "%s" to "%s"', old_id, record.id)

    if len(record.name) > 16:
        record.name = _shorten_ids(record.name)

    if 'accession' in record.annotations and \
       len(record.annotations['accession']) > 16:
        acc = record.annotations['accession']

        record.annotations['accession'] = _shorten_ids(acc)

    # Remove illegal characters from name: otherwise, file cannot be written
    illegal_chars = set('''!"#$%&()*+,:;=>?@[]^`'{|}/ ''')
    for char in record.id:
        if char in illegal_chars:
            record.id = record.id.replace(char, "")
    for char in record.name:
        if char in illegal_chars:
            record.name = record.name.replace(char, "")

def fix_locus_tags(seq_record) -> None: # TODO should be part of secmet
    "Fix CDS feature that don't have a locus_tag, gene name or protein id"
    next_locus_tag = 1

    for feature in seq_record.get_cds_features():
        if feature.get_name() == "no_tag_found":
            logging.critical("fix_locus_tags overwriting tag 'no_tag_found' at %s", feature.location)
            feature.locus_tag = '%s_%05d' % (seq_record.id, next_locus_tag)
            next_locus_tag += 1
        #Fix locus tags, gene names or protein IDs if they contain illegal chars
        illegal_chars = set('''!"#$%&()*+,:; \r\n\t=>?@[]^`'{|}/ ''')
        for attr in ["locus_tag", "gene", "protein_id"]:
            val = getattr(feature, attr)
            if not val or not set(val).intersection(illegal_chars):
                continue
            for char in val:
                if char in illegal_chars:
                    val = val.replace(char, "_")
            logging.critical("%s altered in fix_locus_tags: %s->%s", attr, getattr(feature, attr), val)
            setattr(feature, attr, val)
