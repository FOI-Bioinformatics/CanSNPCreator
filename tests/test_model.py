"""The reference model, the coordinate map, and the parsers that build them.

Each test names the defect it guards against. Before the reference model the
tools held the genome as one concatenated string with a list of cumulative
contig lengths, and converted between the two coordinate systems with running
``correction`` and ``chromosome`` variables at every call site.
"""

import pytest

from cansnpscorer.cansnp_evaluator import annotations_by_contig, locate_snps
from cansnpscorer.helpers.file_utils import (
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
from cansnpscorer.vcf import Snp

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


def test_a_reference_row_longer_than_the_reference_is_refused():
    reference = Reference.from_records([("c1", "ACGTACGTAC")])
    with pytest.raises(ReferenceError, match="once its gaps are removed"):
        AlignmentMap("ACGTACGTACGT", reference)


# -- alignments projected onto reference coordinates -------------------------
#
# x2fa, snippy and RepGenR2 write every row in reference coordinates. Where the
# reference was not aligned the reference row holds a gap at the same column,
# and the row may stop before the reference ends.


def test_a_projected_row_maps_each_column_onto_its_own_position():
    reference = Reference.from_records([("c1", "ACGTA"), ("c2", "CGTAC")])
    mapping = AlignmentMap("AC--ACG-A", reference)
    assert mapping.projected
    assert [mapping.to_reference(column) for column in range(9)] == [
        0, 1, None, None, 4, 5, 6, None, 8]
    assert mapping.locate(5) == (reference.contig("c2"), 1)
    assert mapping.unaligned == 4, "two gaps, one more, and the last base"


def test_a_position_is_found_in_either_layout():
    reference = Reference.from_records([("c1", "ACGTACGTAC")])
    gapped = AlignmentMap("AC--GTACGTAC", reference)
    projected = AlignmentMap("AC-TACGT", reference)
    assert not gapped.projected and projected.projected
    assert gapped.to_column(2) == 4
    assert projected.to_column(2) is None, "that position is not aligned"
    assert projected.to_column(3) == 3
    assert projected.to_column(9) is None, "past the end of the row"
    with pytest.raises(ReferenceError):
        projected.to_column(10)


def test_a_projected_row_of_another_genome_is_refused():
    reference = Reference.from_records([("c1", "ACGTACGTAC")])
    with pytest.raises(ReferenceError, match="c1:2 has G in the alignment"):
        AlignmentMap("AGG-A", reference)


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
    sequence = dataset.sequences["ref"]
    moved = []
    for position in (100, 200, 300):
        others = [base for base in "ACGT" if base != sequence[position - 1]]
        moved.append(Snp("c1", position, "A.1", others[0], others[1]))
    with pytest.raises(ValueError, match="do not fit the reference genome"):
        locate_snps(moved, reference)


def test_snps_outside_the_contig_are_refused(dataset):
    reference = reference_parser(dataset.reference)
    with pytest.raises(ValueError, match="outside the contig"):
        locate_snps([Snp("c1", 10_000, "A.1", "G", "C")], reference)


def test_a_snp_on_an_unknown_contig_is_refused(dataset):
    reference = reference_parser(dataset.reference)
    with pytest.raises(ValueError, match="no contig 'plasmid'"):
        locate_snps([Snp("plasmid", 1, "A.1", "G", "C")], reference)


def test_annotations_are_matched_to_contigs_by_name(dataset):
    reference = reference_parser(dataset.reference)
    matched = annotations_by_contig({"c2": [(1, 2)], "c1": [(3, 4)]}, reference, "gene")
    assert matched == {"c1": [(3, 4)], "c2": [(1, 2)]}

    # A contig the annotations do not mention is empty rather than missing.
    assert annotations_by_contig({"c1": [(3, 4)]}, reference, "gene")["c2"] == []


def test_annotations_whose_names_match_nothing_are_refused(dataset):
    # Matching by file order once attached the genes of a sequence the
    # reference does not hold to another contig, with only a warning.
    reference = reference_parser(dataset.reference)
    with pytest.raises(ValueError, match="none of the 2 sequences"):
        annotations_by_contig({"NC0001": [(1, 2)], "NC0002": [(3, 4)]}, reference, "gene")


def test_annotations_on_a_sequence_the_reference_lacks_are_left_out(dataset, caplog):
    reference = reference_parser(dataset.reference)
    matched = annotations_by_contig({"c1": [(3, 4)], "CP000890.1": [(1, 2)]}, reference, "gene")
    assert matched == {"c1": [(3, 4)], "c2": []}
    assert "CP000890.1" in caplog.text


def test_a_contig_is_named_by_the_first_word_of_its_header(dataset):
    # The dataset's headers read "c1 Testus organismus c1, complete sequence".
    assert reference_parser(dataset.reference).names == ["c1", "c2"]


# -- standard names and annotations ------------------------------------------


def test_two_records_with_one_id_are_refused(tmp_path):
    from cansnpscorer.helpers.file_utils import read_fasta

    path = tmp_path / "twice.fasta"
    path.write_text(">c1 chromosome\nACGT\n>c1 plasmid\nACGT\n")
    with pytest.raises(ValueError, match="names c1 more than once"):
        read_fasta(str(path))


def test_a_prokka_gff3_with_cds_only_gives_genes_and_frames(tmp_path):
    from cansnpscorer.helpers.file_utils import read_gff3

    # Prokka writes CDS features without gene features unless --addgenes,
    # and ends the file with the sequence after a ##FASTA directive.
    path = tmp_path / "prokka.gff"
    path.write_text("##gff-version 3\n"
                    "c1\tProdigal:002006\tCDS\t3\t11\t.\t-\t0\tID=X_00001\n"
                    "c1\tAragorn:001002\ttRNA\t20\t90\t.\t+\t.\tID=X_00002\n"
                    "##FASTA\n>c1\nACGTACGTACGT\n")
    model = read_gff3(str(path))
    assert model.genes == {"c1": [(2, 10, "-")]}
    assert model.cds == {"c1": [(2, 10, "-", 0)]}


@pytest.mark.parametrize("line,message", [
    ("c1\tx\tCDS\t3\t11\t.\t+\t.\tID=a", "without a phase"),
    ("c1\tx\tgene\t11\t3\t.\t+\t.\tID=a", "impossible interval"),
    ("c1 x gene 3 11 . + . ID=a", "9"),
])
def test_a_malformed_gff3_line_is_refused(tmp_path, line, message):
    from cansnpscorer.helpers.file_utils import read_gff3

    path = tmp_path / "bad.gff"
    path.write_text("##gff-version 3\n" + line + "\n")
    with pytest.raises(ValueError, match=message):
        read_gff3(str(path))


def test_genome_names_are_resolved_to_alignment_records():
    from cansnpscorer.helpers.names import resolve_clusters

    records = ["ref", "Fam_Gen_sp_GCF_000000001.1", "strain-7"]
    clusters = {"A": ["GCF_000000001.1", "strain_7", "strain7"], "B": ["missing"]}
    resolved, unmatched, ambiguous = resolve_clusters(clusters, records)
    assert resolved == {"A": ["Fam_Gen_sp_GCF_000000001.1", "strain-7"], "B": []}
    assert unmatched == ["missing"] and ambiguous == []
