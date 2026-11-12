"""The reference model, the coordinate map, and the parsers that build them.

Each test names the defect it guards against. Before the reference model the
tools held the genome as one concatenated string with a list of cumulative
contig lengths, and converted between the two coordinate systems with running
``correction`` and ``chromosome`` variables at every call site.
"""

import pytest

from cansnpscorer.cansnp_evaluator import annotations_by_contig, locate_snps
from cansnpscorer.helpers.file_utils import (
    ParsedSnp,
    fasta_parser,
    read_alignment,
    read_fasta,
    reference_parser,
)
from cansnpscorer.model import (
    AlignmentMap,
    Contig,
    Reference,
    ReferenceError,
    reverse_complement,
    ungap,
)

# -- contigs and references --------------------------------------------------


def write(path, text):
    path.write_text(text)
    return str(path)


def test_a_single_contig_reference_is_usable(tmp_path):
    # The parser recorded a contig's length only when it met the next header,
    # so a one contig reference reported no lengths at all and the evaluator
    # raised IndexError on reference[1] before scoring anything.
    path = write(tmp_path / "ref.fasta", ">chr1\n" + "ACGT" * 100 + "\n")
    reference = reference_parser(path)
    assert len(reference) == 1
    assert reference.length == 400
    assert reference.names == ["chr1"]
    contig, local_position = reference.locate(399)
    assert (contig.name, local_position) == ("chr1", 399)


def test_every_contig_including_the_last_is_measured(tmp_path):
    path = write(
        tmp_path / "ref.fasta",
        ">c1\n" + "A" * 100 + "\n>c2\n" + "C" * 50 + "\n>c3\n" + "G" * 25 + "\n",
    )
    reference = reference_parser(path)
    assert [contig.length for contig in reference] == [100, 50, 25]
    assert reference.length == 175


def test_positions_resolve_on_the_third_contig(tmp_path):
    # The evaluator advanced its correction but never its chromosome_length, so
    # from the third contig on it placed SNPs on the wrong sequence.
    path = write(
        tmp_path / "ref.fasta",
        ">c1\n" + "A" * 100 + "\n>c2\n" + "C" * 50 + "\n>c3\n" + "G" * 25 + "\n",
    )
    reference = reference_parser(path)
    for position, expected in [(0, ("c1", 0)), (99, ("c1", 99)), (100, ("c2", 0)),
                               (149, ("c2", 49)), (150, ("c3", 0)), (174, ("c3", 24))]:
        contig, local_position = reference.locate(position)
        assert (contig.name, local_position) == expected, position


def test_a_position_outside_the_reference_is_refused(tmp_path):
    path = write(tmp_path / "ref.fasta", ">c1\nACGT\n")
    reference = reference_parser(path)
    with pytest.raises(ReferenceError, match="outside the reference"):
        reference.locate(4)


def test_a_repeated_contig_name_is_refused():
    with pytest.raises(ReferenceError, match="more than once"):
        Reference.from_records([("c1", "ACGT"), ("c1", "TTTT")])


def test_the_window_around_a_position_is_clamped_to_the_contig():
    # Slicing with a negative start wrapped to the end of the sequence, which
    # gave a window that was short, in the wrong place, or empty.
    contig = Contig(name="c1", sequence="ACGT" * 50, start=0)
    assert contig.window(10, 50) == contig.sequence[0:60]
    assert len(contig.window(0, 50)) == 50
    assert len(contig.window(199, 50)) == 51


# -- the alignment to reference map ------------------------------------------


def test_columns_map_straight_through_without_gaps():
    reference = Reference.from_records([("c1", "ACGTACGTAC")])
    mapping = AlignmentMap("ACGTACGTAC", reference)
    assert [mapping.to_reference(column) for column in range(10)] == list(range(10))
    assert mapping.locate(9) == (reference.contig("c1"), 10)


def test_a_gap_in_the_reference_row_shifts_the_columns_after_it():
    # The extractor reported the alignment column as if it were a reference
    # position, so every gap in the reference row moved all later SNPs.
    reference = Reference.from_records([("c1", "ACGTACGTAA")])
    mapping = AlignmentMap("AC-GTACGTAA", reference)
    assert mapping.to_reference(1) == 1, "before the gap, column and position agree"
    assert mapping.to_reference(2) is None, "the reference is gapped here"
    assert mapping.to_reference(3) == 2, "after the gap they differ by one"
    assert mapping.to_reference(10) == 9
    contig, position = mapping.locate(10)
    assert (contig.name, position) == ("c1", 10), "reported 1-based"


def test_a_reference_row_of_the_wrong_length_is_refused():
    reference = Reference.from_records([("c1", "ACGTACGTAC")])
    with pytest.raises(ReferenceError, match="once its gaps are removed"):
        AlignmentMap("ACGTACGT", reference)


def test_an_alignment_of_another_genome_is_refused():
    reference = Reference.from_records([("c1", "ACGTACGTAC")])
    with pytest.raises(ReferenceError, match="c1:5 has T in the alignment and A in the reference"):
        AlignmentMap("ACGTTCGTAC", reference)


def test_ungap_and_reverse_complement():
    assert ungap("AC-GT.AC~GT") == "ACGTACGT"
    assert reverse_complement("ACGTN") == "NACGT"


# -- the FASTA parsers -------------------------------------------------------


def test_a_repeated_sequence_name_is_refused(tmp_path):
    # The dict built by the parser silently kept only the last of them.
    path = write(tmp_path / "a.fasta", ">ref\nACGT\n>g1\nACGA\n>g1\nTCGT\n")
    with pytest.raises(ValueError, match="names g1 more than once"):
        read_fasta(path)


def test_an_alignment_that_is_not_aligned_is_refused(tmp_path):
    # Indexing each genome by column raised a bare IndexError part way through.
    path = write(tmp_path / "a.fasta", ">ref\nACGTACGT\n>g1\nACGT\n")
    with pytest.raises(ValueError, match="not aligned"):
        read_alignment(path)


def test_an_alignment_needs_more_than_the_reference(tmp_path):
    path = write(tmp_path / "a.fasta", ">ref\nACGTACGT\n")
    with pytest.raises(ValueError, match="at least one other genome"):
        read_alignment(path)


def test_soft_masked_sequence_is_upper_cased(tmp_path):
    # Lower case bases were counted as unusable, so SNPs in soft masked regions
    # disappeared and the GC content of such a region read as zero.
    path = write(tmp_path / "a.fasta", ">c1\nacgtACGT\n")
    assert read_fasta(path) == [("c1", "ACGTACGT")]


def test_sequence_before_the_first_header_is_refused(tmp_path):
    path = write(tmp_path / "a.fasta", "ACGT\n>c1\nACGT\n")
    with pytest.raises(ValueError, match="before the first header"):
        read_fasta(path)


def test_fasta_parser_keeps_its_three_output_formats(tmp_path):
    path = write(tmp_path / "a.fasta", ">c1\nACGT\n>c2\nTTTT\n")
    assert fasta_parser(path, "dict") == {"c1": "ACGT", "c2": "TTTT"}
    assert fasta_parser(path, "named_list") == [("c1", "ACGT"), ("c2", "TTTT")]
    assert fasta_parser(path, "list") == ["ACGT", "TTTT"]


# -- placing SNPs against the reference --------------------------------------


def test_snps_are_checked_against_the_reference_base(dataset, parsed_snps):
    reference = reference_parser(dataset.reference)
    located = locate_snps(parsed_snps, reference)
    assert len(located) == len(parsed_snps)
    for contig, local_position, _node, alleles in located:
        assert contig.base(local_position) in alleles


def test_snps_from_another_reference_are_refused(dataset):
    # The direct test for a coordinate mismatch: positions produced against a
    # different reference, or against alignment columns, disagree with the base
    # the reference carries.
    reference = reference_parser(dataset.reference)
    moved = [ParsedSnp(node="A.1", alleles=("G", "C"), contig="c1", position=position + 7)
             for position in (100, 200, 300)]
    with pytest.raises(ValueError, match="disagree with the reference genome"):
        locate_snps(moved, reference)


def test_snps_outside_the_reference_are_refused(dataset):
    reference = reference_parser(dataset.reference)
    beyond = ParsedSnp(node="A.1", alleles=("G", "C"), contig="c1", position=10_000)
    with pytest.raises(ValueError, match="outside the reference"):
        locate_snps([beyond], reference)

    legacy = ParsedSnp(node="A.1", alleles=("G", "C"), global_position=reference.length)
    with pytest.raises(ValueError, match="outside the reference"):
        locate_snps([legacy], reference)


def test_a_snp_on_an_unknown_contig_is_refused(dataset):
    reference = reference_parser(dataset.reference)
    unknown = ParsedSnp(node="A.1", alleles=("G", "C"), contig="plasmid", position=1)
    with pytest.raises(ValueError, match="no contig 'plasmid'"):
        locate_snps([unknown], reference)


def test_positions_from_the_older_format_are_placed_on_their_contig(dataset, cansnps):
    # The older format counts over the concatenated genome from zero.
    reference = reference_parser(dataset.reference)
    on_second = next(cansnp for cansnp in cansnps if cansnp.column >= 700)
    legacy = ParsedSnp(node=on_second.node, alleles=on_second.bases,
                       global_position=on_second.column)
    contig, local_position, _node, _alleles = locate_snps([legacy], reference)[0]
    assert (contig.name, local_position) == ("c2", on_second.column - 700)


def test_annotations_must_describe_the_same_contigs_as_the_reference(dataset):
    reference = reference_parser(dataset.reference)
    assert sorted(annotations_by_contig({1: [], 2: []}, reference, "GFF")) == ["c1", "c2"]
    with pytest.raises(ValueError, match="describes 1 contigs but the reference has 2"):
        annotations_by_contig({1: []}, reference, "GFF")
