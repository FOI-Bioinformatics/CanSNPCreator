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
    for contig, local_position, _node, alleles, _nearest in located:
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
    assert model.genes == {"c1": [(2, 10)]}
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



# -- found in review ---------------------------------------------------------


def test_record_ids_with_a_fasta_suffix_match_the_genome_map():
    # Parsnp names alignment records after their files, as g1.fna.
    from cansnpscorer.helpers.names import resolve_names

    matched, unmatched, _ambiguous = resolve_names(["g1.fna", "g2"], ["g1.fna", "g2.fna"])
    assert matched == {"g1.fna": "g1.fna", "g2": "g2.fna"} and unmatched == []


def test_an_exact_name_wins_and_a_collision_is_reported():
    from cansnpscorer.helpers.names import resolve_names

    matched, _unmatched, ambiguous = resolve_names(["A_1", "A1"], ["A_1", "A-1"])
    assert matched == {"A_1": "A_1"}
    assert ambiguous == ["A1"], "A1 could be either record"


def test_a_snp_in_a_gene_holding_smaller_features_is_in_a_gene(tmp_path):
    from cansnpscorer.helpers.file_utils import read_gff3
    from cansnpscorer.helpers.scoring_methods import SnpContext, gene_score
    from cansnpscorer.model import Contig

    path = tmp_path / "nested.gff"
    path.write_text("##gff-version 3\n" + "".join(
        f"c1\tx\tgene\t{start}\t{end}\t.\t+\t.\tID=g{start}\n"
        for start, end in [(1, 1000), (11, 20), (31, 40), (51, 60), (71, 80)]))
    genes = read_gff3(str(path)).genes["c1"]
    assert genes == [(0, 999)]
    snp = SnpContext(Contig("c1", "A" * 1000, 0), 500, "A", ("A", "G"), genes=genes)
    assert gene_score(snp, 0) == 1.0


def test_a_codon_holding_an_unknown_base_has_no_amino_acid():
    from cansnpscorer.helpers.scoring_methods import SnpContext, mutation_score
    from cansnpscorer.model import Contig

    snp = SnpContext(Contig("c1", "ATGNAATAA", 0), 4, "A", ("C", "A"), cds=[(0, 8, "+", 0)])
    assert mutation_score(snp) == 0.0


def test_a_later_cds_is_read_when_the_first_cannot_give_a_codon():
    from cansnpscorer.helpers.scoring_methods import SnpContext, mutation_score
    from cansnpscorer.model import Contig

    # In the first CDS the SNP's codon would run past the CDS's end; the
    # second reads it whole.
    snp = SnpContext(Contig("c1", "TTATTGTAA", 0), 2, "A", ("A", "G"),
                     cds=[(0, 3, "+", 2), (0, 8, "+", 0)])
    assert mutation_score(snp) == 1.0, "TTA and TTG are both Leu"


@pytest.mark.parametrize("line", [
    "##INFO=<Number=1,Type=String,Description=\"no ID\">",
    "##INFO=<ID=X,Number=1,Type=String,Description=\"unterminated>",
    "##contig=<ID=c1,length=many>",
])
def test_a_malformed_vcf_header_line_is_refused_cleanly(tmp_path, line):
    from cansnpscorer import vcf

    path = tmp_path / "bad.vcf"
    path.write_text("##fileformat=VCFv4.3\n" + line + "\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    with pytest.raises(vcf.VcfError, match="malformed header line"):
        vcf.read(str(path))


def test_an_empty_bed_region_is_skipped(tmp_path):
    from cansnpscorer.helpers.file_utils import read_intervals

    path = tmp_path / "r.bed"
    path.write_text("c1\t5\t5\nc1\t10\t20\n")
    assert read_intervals(str(path)) == {"c1": [(10, 19)]}


def test_a_root_given_as_its_own_parent_is_the_root(tmp_path):
    from cansnpscorer.helpers.file_utils import read_hierarchy

    path = tmp_path / "nodes.tsv"
    path.write_text("child\tparent\n1\t1\n2\t1\n")
    assert read_hierarchy(str(path)).root == "1"


def test_masks_over_the_concatenated_reference_are_placed_on_contigs():
    from cansnpscorer.cansnp_evaluator import regions_over_contigs

    reference = Reference.from_records([("c1", "A" * 700), ("c2", "A" * 300)])
    # Gubbins reports alignment columns whatever the sequence is called; one
    # region here spans the end of c1 and the start of c2.
    placed = regions_over_contigs({"SEQUENCE": [(10, 19), (690, 709)]}, reference)
    assert placed == {"c1": [(10, 19), (690, 699)], "c2": [(0, 9)]}
    with pytest.raises(ValueError, match="past the end of the reference"):
        regions_over_contigs({"SEQUENCE": [(995, 1005)]}, reference)


# -- annotations from Prokka, Bakta and NCBI --------------------------------


def two_contig_reference():
    return Reference.from_records([("NZ_A.1", "ATGAAACCCGGGTAA" * 4), ("NZ_B.1", "GGGCCCAAATTT" * 3)])


def test_contigs_renamed_by_prokka_are_matched_by_their_sequence(tmp_path):
    # prokka --compliant renames contigs gnl|centre|PREFIX_n, but writes the
    # sequences after ##FASTA, by which they are matched.
    from cansnpscorer.cansnp_evaluator import annotations_by_contig
    from cansnpscorer.helpers.file_utils import read_gff3

    reference = two_contig_reference()
    path = tmp_path / "prokka.gff"
    path.write_text("##gff-version 3\n"
                    "gnl|X|PFX_1\tProdigal:002006\tCDS\t1\t15\t.\t+\t0\tID=PFX_00001\n"
                    "##FASTA\n>gnl|X|PFX_1\n" + "atgaaacccgggtaa" * 4 + "\n"
                    ">gnl|X|PFX_2\n" + "GGGCCCAAATTT" * 3 + "\n")
    model = read_gff3(str(path))
    assert set(model.sequences) == {"gnl|X|PFX_1", "gnl|X|PFX_2"}
    genes = annotations_by_contig(model.genes, reference, "gene", sequences=model.sequences)
    assert genes == {"NZ_A.1": [(0, 14)], "NZ_B.1": []}


def test_a_sequence_that_differs_is_not_matched(tmp_path):
    from cansnpscorer.cansnp_evaluator import annotations_by_contig

    reference = two_contig_reference()
    with pytest.raises(ValueError, match="--contig-aliases"):
        annotations_by_contig({"other": [(0, 1)]}, reference, "gene",
                              sequences={"other": "ATGAAACCCGGGTAG" * 4})


def test_an_ncbi_assembly_report_pairs_genbank_and_refseq_names(tmp_path):
    from cansnpscorer.cansnp_evaluator import annotations_by_contig
    from cansnpscorer.helpers.file_utils import read_contig_aliases

    report = tmp_path / "assembly_report.txt"
    report.write_text(
        "# Assembly name:  ASM15779v1\n"
        "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\t"
        "GenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\t"
        "UCSC-style-name\n"
        "supercont1.1\tunplaced-scaffold\tna\tna\tA.1\t=\tNZ_A.1\tPrimary Assembly\t60\tna\n"
        "supercont1.2\tunplaced-scaffold\tna\tna\tB.1\t=\tNZ_B.1\tPrimary Assembly\t36\tna\n")
    aliases = read_contig_aliases(str(report))
    assert {"supercont1.1", "A.1", "NZ_A.1"} in aliases
    genes = annotations_by_contig({"A.1": [(0, 5)], "B.1": [(1, 2)]}, two_contig_reference(),
                                  "gene", aliases)
    assert genes == {"NZ_A.1": [(0, 5)], "NZ_B.1": [(1, 2)]}


def test_a_ucsc_chromalias_table_is_read(tmp_path):
    from cansnpscorer.helpers.file_utils import read_contig_aliases

    table = tmp_path / "chromAlias.txt"
    table.write_text("# refseq\tgenbank\tucsc\nNZ_A.1\tA.1\tchr1\n")
    assert read_contig_aliases(str(table)) == [{"NZ_A.1", "A.1", "chr1"}]


def test_two_names_for_one_contig_are_refused():
    from cansnpscorer.cansnp_evaluator import annotations_by_contig

    with pytest.raises(ValueError, match="describes the contig NZ_A.1 twice"):
        annotations_by_contig({"NZ_A.1": [(0, 1)], "A.1": [(2, 3)]}, two_contig_reference(),
                              "gene", [{"A.1", "NZ_A.1"}])


def test_pseudogenes_and_other_genetic_codes_give_no_frame(tmp_path, caplog):
    from cansnpscorer.helpers.file_utils import read_gff3

    path = tmp_path / "ncbi.gff"
    path.write_text(
        "##gff-version 3\n"
        "NZ_A.1\tRefSeq\tCDS\t1\t15\t.\t+\t0\tID=cds-1;pseudo=true;transl_table=11\n"
        "NZ_A.1\tRefSeq\tCDS\t16\t30\t.\t+\t0\tID=cds-2;transl_table=4\n"
        "NZ_A.1\tRefSeq\tCDS\t31\t45\t.\t+\t0\tID=cds-3;product=a%3Bb;transl_table=11\n")
    model = read_gff3(str(path))
    assert model.cds == {"NZ_A.1": [(30, 44, "+", 0)]}
    assert model.genes == {"NZ_A.1": [(0, 44)]}, "all three still place a SNP in a gene"
    assert "genetic code 4" in caplog.text


def test_gzipped_fasta_and_gff3_are_read(tmp_path):
    import gzip

    from cansnpscorer.helpers.file_utils import read_fasta, read_gff3

    fasta = tmp_path / "ref.fna.gz"
    with gzip.open(fasta, "wt") as handle:
        handle.write(">NZ_A.1 description\nACGT\n")
    gff = tmp_path / "genes.gff.gz"
    with gzip.open(gff, "wt") as handle:
        handle.write("##gff-version 3\nNZ_A.1\tx\tCDS\t1\t3\t.\t+\t0\tID=a\n")
    assert read_fasta(str(fasta)) == [("NZ_A.1", "ACGT")]
    assert read_gff3(str(gff)).cds == {"NZ_A.1": [(0, 2, "+", 0)]}


# -- the variation around a SNP ----------------------------------------------


def test_variation_is_counted_among_the_chosen_genomes_only():
    from cansnpscorer.helpers.data_utils import variable_columns

    alignment = {
        "out": "AAAAAAAA",
        "x": "ACAAAA-A",
        "y": "AAGAAAAN",
        "z": "AAAATAAA",
    }
    # Column 1 and 2 differ among x, y and z, a third base included; column 4
    # differs in one genome; column 6 only by a gap and 7 by an N, which are
    # no evidence of a different base.
    assert variable_columns(alignment, ["x", "y", "z"]) == [1, 2, 4]
    # A genome left out, as the outgroup is, does not count: without z,
    # column 4 no longer differs.
    assert variable_columns(alignment, ["out", "x", "y"]) == [1, 2]
    assert variable_columns(alignment, ["x"]) == []


def test_the_nearest_variation_is_sought_on_the_same_contig():
    from cansnpscorer.helpers.data_utils import nearest_distances

    variable = {"c1": [10, 15, 40], "c2": [3]}
    assert nearest_distances([("c1", 15), ("c1", 12), ("c1", 100), ("c2", 3), ("c3", 5)],
                             variable) == [5, 2, 60, None, None]
