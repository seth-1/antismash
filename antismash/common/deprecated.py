# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
This file will be removed as soon as the new abstraction layer is completed.
"""

import logging
import os
import re
import sys

import Bio
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation # for others importing
from Bio.SeqRecord import SeqRecord
from helperlibs.bio import seqio

from antismash.common import gff_parser
from antismash.common.secmet import Record, CDSFeature

# temporary code skip logging # TODO
import inspect
import linecache

def CODE_SKIP_WARNING():
    prev = inspect.currentframe().f_back
    logging.critical("skipping code:" + prev.f_code.co_name +"():" \
            + linecache.getline(prev.f_code.co_filename, prev.f_lineno + 1).replace('%', '%%'))
# end temp


def get_all_features_of_type(seq_record, types):
    "Return all features of the specified types for a seq_record"
    logging.critical("utils.get_all_features_of_type() called")
    raise RuntimeError("get_all_features_of_type(record, types) called, did you mean record.get_*()")

def get_cds_features_within_clusters(seq_record):
    cds_features = []
    for cluster in seq_record.get_clusters():
        cds_features.extend(cluster.cds_children)
    return cds_features

def get_withincluster_cds_features(seq_record):
    logging.critical("get_withincluster_cds_features() deprecated, use get_cds_features_within_clusters()")
    return get_cds_features_within_clusters(seq_record)

def parse_input_sequence(filename, options):
    "Parse the input sequences from given filename"
    logging.info('Parsing input sequence %r', filename)

    sequences = []
    if not os.path.exists(filename):
        msg = "Sequence file not found: %r" % filename
        logging.error(msg)
        raise ValueError(msg)

    try:
        record_list = list(seqio.parse(filename))
        if not record_list:
            logging.error('No sequence in file %r', filename)
        sequences.extend([rec for rec in record_list if len(rec.seq) > options.minlength or \
            ('contig' in rec.annotations or 'wgs_scafld' in rec.annotations or \
            'wgs' in rec.annotations)])
    except (ValueError, AssertionError) as err:
        logging.error('Parsing %r failed: %s', filename, err)
        raise
    except Exception as err:
        logging.error('Parsing %r failed with unhandled exception: %s',
                      filename, err)
        raise
    return [Record.from_biopython(sequence) for sequence in sequences]

def pre_process_sequences(sequences, options, genefinding):
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
        if not sequence.seq or (
                options.input_type == 'nucl' and \
                not str(sequence.seq).replace("N", "") and \
                cdsfeatures_with_translations < 0.8 * len(cdsfeatures)):
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

    #Handle WGS master or supercontig entries
    check_for_wgs_scaffolds(sequences)

    #Now remove small contigs < minimum length again
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
                sequence.skip = "skipping all but first {0} meaningful records (--limit {0}) ".format(options.limit)

    if warned:
        options.triggered_limit = True

    #Check if no duplicate locus tags / gene IDs are found
    check_duplicate_gene_ids(sequences)

    # Check GFF suitability
    if options.genefinding_gff3:
        single_entry = gff_parser.check_gff_suitability(options, sequences)

    # Ensure all records have valid names
    for seq_record in sequences:
        fix_record_name_id(seq_record, {seq.id for seq in sequences}, options)

    for sequence in sequences:
        if sequence.skip:
            continue
        #Fix sequence name (will be ID) if it contains illegal chars
        illegal_chars = '''!"#$%&()*+,:; \r\n\t=>?@[]^`'{|}/ '''
        for char in sequence.name:
            if char in illegal_chars:
                sequence.name = sequence.name.replace(char, "_")
        #Iterate through sequence objects
        if len(sequence.get_cds_features()) < 1:
            if options.genefinding_gff3:
                logging.info("No CDS features found in record %r but GFF3 file provided, running GFF parser.", sequence.id)
                gff_parser.run(sequence, single_entry, options)
                check_duplicate_gene_ids(sequences)
            elif options.genefinding_tool != "none":
                logging.info("No CDS features found in record %r, running gene finding.", sequence.id)
                genefinding.run_on_record(sequence, options)
            if len(sequence.get_cds_features()) < 1:
                logging.info("No genes found, skipping record")
                sequence.skip = "No genes found"
                continue
        #Fix locus tags
        fix_locus_tags(sequence)

    #Make sure that all CDS entries in all seq_records have translation tags, otherwise add them
    add_translations(sequences)

    #Make sure that all seq_records have a sequence
    add_seq_record_seq(sequences)

    if len(sequences) > 1 and (options.start != -1 or options.end != -1):
        options.start = -1
        options.end = -1
        logging.info("Discarding --start and --end options, as multiple entries are used.")

    for i, sequence in enumerate(sequences):
        if options.start > 1:
            if options.start > len(sequence):
                logging.error('Specified analysis start point is at %r, which is larger ' \
                              'than record size %r', options.start, len(sequence))
                sys.exit(1)
            sequence = sequence[options.start-1:]
            # new sequence is shorter, so fix the end calculation
            options.end -= options.start
            sequences[i] = sequence

        if options.end > 0:
            if options.end > len(sequence):
                logging.error('Specified analysis end point is at %r, which is larger ' \
                              'than record size %r', options.end, len(sequence))
                sys.exit(1)
            sequence = sequence[:options.end]
            sequences[i] = sequence

        # Some programs write gaps as - or X or x, translate requires N
        for char in ['-', 'X', 'x']:
            sequence.seq = Seq(str(sequence.seq).replace(char, 'N'),
                               alphabet=sequence.seq.alphabet)



    #Fix sequence record IDs to be unique
    ids_used = []
    for sequence in sequences:
        seq_id = sequence.id
        if seq_id not in ids_used:
            ids_used.append(seq_id)
            continue
        prefix = seq_id
        if len(prefix) > 11:
            prefix = prefix[:11]


        suffix = 0
        #Make sure the length of the ID does not exceed 16
        if len(seq_id) <= 11:
            while "%s_%i" % (seq_id, suffix) in ids_used:
                suffix += 1
            sequence.id = "%s_%i" % (seq_id, suffix)
        else:
            while "%s_%i" % (seq_id[:-4], suffix) in ids_used:
                suffix += 1
            sequence.id = "%s_%i" % (seq_id[:-4], suffix)
        ids_used.append(sequence.id)
    return sequences

def generate_unique_id(prefix, existing_ids, start=0, max_length=-1):
    """ Generate a identifier of the form prefix_num, e.g. seq_15.

        Args:
            prefix: The text portion of the name.
            existing_ids: The current identifiers to avoid collision with.
            start: An integer to start counting at (default: 0)
            max_length: The maximum length allowed for the identifier,
                        values less than 1 are considerd to be no limit.

        Returns:
            A tuple of the identifier generated and the value of the counter
                at the time the identifier was generated, e.g. ("seq_15", 15)

    """
    counter = int(start)
    existing_ids = set(existing_ids)

    format_string = "{}_{}".format(prefix, counter)
    name = format_string % counter
    while name in existing_ids:
        counter += 1
        name = format_string % counter
    if max_length > 0 and len(name) > max_length:
        raise RuntimeError("Could not generate unique id for %s after %d iterations" % (prefix, counter - start))
    return name, counter

def is_nucl_seq(sequence):
    other = str(sequence).lower()
    for char in "acgtn":
        other = other.replace(char, "")
    return len(other) < 0.2 * len(sequence)

def generate_nucl_seq_record(sequences):
    "Generate single nucleotide seq_record from supplied sequences"
    if not sequences:
        raise ValueError("Cannot generate nucleotide records of empty input")
    record = Record(Seq(""), id="Protein_Input", name="ProteinInput",
                   description="antiSMASH protein input")
    position = 0
    cdsnames = set()
    for sequence in sequences:
        startpos = position
        endpos = position + len(sequence) * 3
        position += len(sequence) * 3 + 1000
        location = FeatureLocation(startpos, endpos, strand=1)
        name = sequence.id[:15].replace(" ", "_")
        if name in cdsnames:
            name, _ = generate_unique_id(name[:8], cdsnames)
        cdsnames.add(name)
        translation = str(sequence.seq).replace('.', 'X')
        cdsfeature = CDSFeature(location, translation, product=name, locus_tag=name)
        record.add_cds_feature(cdsfeature)
    return record

def check_duplicate_gene_ids(sequences):
    "Fix duplicate locus tags so that they are different"
    no_tag = "no_tag_found"
    high_water_mark = 0
    all_ids = set()
    for sequence in sequences:
        for cdsfeature in sequence.get_cds_features():
            name = get_gene_id(cdsfeature)
            if not name:
                name = no_tag
            if name == no_tag:
                name, high_water_mark = generate_unique_id(name[:8], all_ids,
                                                start=high_water_mark + 1)
            elif name in all_ids:
                name, _ = generate_unique_id(name[:8], all_ids, start=1)
            cdsfeature.product = name
            cdsfeature.locus_tag = name
            all_ids.add(name)


def fix_locus_tags(seq_record):
    "Fix CDS feature that don't have a locus_tag, gene name or protein id"
    next_locus_tag = 1

    for feature in seq_record.get_cds_features():
        if get_gene_id(feature) == "no_tag_found":
            feature.locus_tag = 'AUTOORF_%05d' % next_locus_tag
            next_locus_tag += 1
        #Fix locus tags, gene names or protein IDs if they contain illegal chars
        illegal_chars = '''!"#$%&()*+,:; \r\n\t=>?@[]^`'{|}/ '''
        for attr in ["locus_tag", "gene", "protein_id"]:
            val = getattr(feature, attr)
            if not val:
                continue
            for char in val:
                if char in illegal_chars:
                    val = val.replace(char, "_")
            setattr(feature, attr, val)

def add_translations(seq_records):
    "Add a translation qualifier to all CDS features"
    for seq_record in seq_records:
        if seq_record.skip:
            continue
        logging.debug("Adding translations to record: %s", seq_record.id)
        cdsfeatures = seq_record.get_cds_features()
        for cdsfeature in cdsfeatures:
            if cdsfeature.translation:
                continue
            if not seq_record.seq:
                logging.error('No amino acid sequence in input entry for CDS %r, ' \
                        'and no nucleotide sequence provided to translate it from.', cdsfeature.id)
                raise ValueError("Missing sequence info for CDS %r" % cdsfeature.id)
            try:
                translation = str(get_aa_translation(seq_record, cdsfeature))
            except Bio.Data.CodonTable.TranslationError as err:
                logging.error('Getting amino acid sequences from %s, CDS %r failed: %s',
                        seq_record.name, cdsfeature.id, err)
                raise
            cdsfeature.translation = translation

def get_gene_id(feature):
    "Get the gene ID from locus_tag, gene name or protein id, in that order"
    gene_id = "no_tag_found"
    for label in ['locus_tag', 'gene', 'protein_id']:
        if hasattr(feature, label):
            value = getattr(feature, label)
            if value:
                gene_id = value
                break
    assert isinstance(gene_id, str), type(gene_id)
    return gene_id

def add_seq_record_seq(seq_records):
    for seq_record in seq_records:
        if not seq_record.seq:
            cds_features = seq_record.get_cds_features()
            start_max = max([cds.location.start for cds in cds_features])
            end_max = max([cds.location.end for cds in cds_features])
            seq_record.seq = Seq(max([start_max, end_max]) * "n")

def get_aa_translation(seq_record, feature):
    """Obtain content for translation qualifier for specific CDS feature in sequence record"""
    fasta_seq = feature.extract(seq_record.seq).ungap('-').translate(to_stop=True)
    if not fasta_seq:
        logging.debug("Retranslating %s with stop codons", feature.id)
        fasta_seq = feature.extract(seq_record.seq).ungap('-').translate()
    if "*" in str(fasta_seq):
        fasta_seq = Seq(str(fasta_seq).replace("*", "X"), Bio.Alphabet.generic_protein)
    if "-" in str(fasta_seq):
        fasta_seq = Seq(str(fasta_seq).replace("-", ""), Bio.Alphabet.generic_protein)

    return fasta_seq


def check_for_wgs_scaffolds(seq_records):
    for seq_record in seq_records:
        #Check if seq_record is a WGS master record or a supercontig record
        if 'wgs_scafld' in seq_record.annotations \
                or 'wgs' in seq_record.annotations \
                or 'contig' in seq_record.annotations:
            raise RuntimeError("Incomplete whole genome shotgun records not supported")

def get_feature_dict(seq_record):
    """Get a dictionary mapping features to their IDs"""
    features = seq_record.get_cds_features()
    feature_by_id = {}
    for feature in features:
        gene_id = get_gene_id(feature)
        feature_by_id[gene_id] = feature
    return feature_by_id


def get_multifasta(seq_record):
    """Extract multi-protein FASTA from all CDS features in sequence record"""
    features = seq_record.get_cds_features()
    all_fastas = []
    for feature in features:
        gene_id = get_gene_id(feature)
        fasta_seq = feature.translation
        if "-" in str(fasta_seq):
            fasta_seq = Seq(str(fasta_seq).replace("-", ""), Bio.Alphabet.generic_protein)

        # Never write empty fasta entries
        if not fasta_seq:
            logging.error("No translation for CDS %s", gene_id)
            raise ValueError("No translation for CDS %s" % gene_id)

        all_fastas.append(">%s\n%s" % (gene_id, fasta_seq))
    full_fasta = "\n".join(all_fastas)
    return full_fasta

def get_cluster_type(cluster):
    "Get product type of a gene cluster"
    logging.critical("utils.get_cluster_type() called")
    raise RuntimeError("utils.get_cluster_type(cluster) called, did you mean cluster.product?")

def get_cluster_cds_features(cluster, seq_record):
    logging.critical("utils.get_cluster_cds_features() called")
    raise RuntimeError("utils.get_cluster_cds_features(cluster) called, did you mean cluster.cds_children?")

def strip_record(seq_record):
    """ Discard antismash specific features and feature qualifiers """
    seq_record.clear_clusters()
    seq_record.clear_cluster_borders()
    seq_record.clear_cds_motifs()
    seq_record.clear_antismash_domains()

    # clean up antiSMASH annotations in CDS features
    for feature in seq_record.get_cds_features():
        feature.sec_met = None

def sort_features(seq_record):
    "Sort features in a seq_record by their position"
    logging.critical("utils.sort_features() called")
    raise RuntimeError("utils.sort_features(seq_record) called, did you mean sorted(seq_record.get_all_features())?")

def fix_record_name_id(seq_record, all_record_ids, options):
    "Fix a seq record's name and id to be <= 16 characters, the GenBank limit; if record name is too long, add c000X prefix"

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
            contig_no = seq_record.record_index

        return "c{ctg:05d}_{origid}..".format(ctg=contig_no, origid=idstring[:7])

    if seq_record.id == "unknown.1":
        seq_record.id = "unk_seq_{ctg:05d}".format(ctg=seq_record.record_index)
        logging.warning('Invalid sequence id "unknown.1", replaced by %s', seq_record.id)

    if seq_record.name == "unknown":
        seq_record.name = "unk_seq_{ctg:05d}".format(ctg=options.record_index)
        logging.warning('Invalid sequence name "unknown", replaced by %s', seq_record.name)

    if len(seq_record.id) > 16:
        oldid = seq_record.id

        #Check if it is a RefSeq accession number like NZ_AMZN01000079.1 that is just too long because of the version number behind the dot
        if (seq_record.id[-2] == "." and
                seq_record.id.count(".") == 1 and
                len(seq_record.id.partition(".")[0]) <= 16 and
                seq_record.id.partition(".")[0] not in all_record_ids):
            seq_record.id = seq_record.id.partition(".")[0]
            all_record_ids.add(seq_record.id)
        else: #Check if the ID suggested by _shorten_ids is unique
            if _shorten_ids(oldid) not in all_record_ids:
                name = _shorten_ids(oldid)
            else:
                name, _ = generate_unique_id(seq_record.id[:12], all_record_ids, max_length=16)
            seq_record.id = name
            all_record_ids.add(name)

        logging.warning('Fasta header too long: renamed "%s" to "%s"', oldid, seq_record.id)

    if len(seq_record.name) > 16:

        seq_record.name = _shorten_ids(seq_record.name)

    if 'accession' in seq_record.annotations and \
       len(seq_record.annotations['accession']) > 16:
        acc = seq_record.annotations['accession']

        seq_record.annotations['accession'] = _shorten_ids(acc)

    # Remove illegal characters from name: otherwise, file cannot be written
    illegal_chars = '''!"#$%&()*+,:;=>?@[]^`'{|}/ '''
    for char in seq_record.id:
        if char in illegal_chars:
            seq_record.id = seq_record.id.replace(char, "")
    for char in seq_record.name:
        if char in illegal_chars:
            seq_record.name = seq_record.name.replace(char, "")
