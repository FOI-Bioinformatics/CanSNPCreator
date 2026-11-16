"""One readable report of an audit.

The audit writes its findings as a table and the typing of the genomes as
five more. Those are for further processing. This module puts what a reader
asks first into one Markdown file: what was run, what the scheme holds, what
was found, and how the scheme typed the genomes. It adds no number of its own;
every figure is taken from the findings or from the typing report, so the
tables remain the record.
"""

from collections.abc import Sequence

from cansnpscorer.findings import PROBLEM, RECOMMENDATION, Finding
from cansnpscorer.scheme import Scheme
from cansnpscorer.typing_report import (
    GENOME_COLUMNS,
    ISOLATE,
    NA,
    RANK_COLUMNS,
    SNP_COLUMNS,
    SUMMARY_COLUMNS,
    TypingReport,
)

#: At most this many findings are written out under one category.
LISTED = 10
#: At most this many rows are written to a table of the worst nodes, canSNPs
#: or genomes.
WORST = 15


def cell(text: object) -> str:
    """A value as the cell of a Markdown table."""
    return str(text).replace("|", "\\|").replace("\n", " ")


def table(header: Sequence[str], rows: Sequence[Sequence[object]]) -> list[str]:
    lines = ["| " + " | ".join(cell(name) for name in header) + " |",
             "|" + "|".join(" --- " for _name in header) + "|"]
    lines += ["| " + " | ".join(cell(value) for value in row) + " |" for row in rows]
    return [*lines, ""]


def _scheme_lines(scheme: Scheme) -> list[str]:
    per_node: dict[str, int] = {}
    for snp in scheme.checked_snps:
        per_node[snp.node] = per_node.get(snp.node, 0) + 1
    sizes = sorted(per_node.values())
    rows = [
        ("nodes", len(scheme.nodes)),
        ("nodes with a canSNP", len(per_node)),
        ("canSNPs called by CanSNPer2, one per node", len(scheme.snps)),
        ("further ranked canSNPs", sum(1 for snp in scheme.snp_set if snp.rank > 1)),
        ("canSNPs per node, fewest to most", f"{sizes[0]}-{sizes[-1]}" if sizes else NA),
        ("reference genomes", len(scheme.references)),
        ("genomes listed by the scheme", len(scheme.genomes) or NA),
    ]
    return ["## The scheme", "", *table(("item", "value"), rows)]


def _finding_lines(findings: Sequence[Finding]) -> list[str]:
    lines = ["## Findings", ""]
    if not findings:
        return [*lines, "Nothing was found that asks for a change.", ""]
    counts: dict[tuple[str, str, str], int] = {}
    for finding in findings:
        key = (finding.severity, finding.category, finding.action)
        counts[key] = counts.get(key, 0) + 1
    order = {PROBLEM: 0, RECOMMENDATION: 1}
    ranked = sorted(counts.items(), key=lambda item: (order[item[0][0]], -item[1], item[0][1]))
    lines += table(("severity", "category", "action", "findings"),
                   [(*key, found) for key, found in ranked])
    for severity, title in ((PROBLEM, "Problems"), (RECOMMENDATION, "Recommendations")):
        chosen = [finding for finding in findings if finding.severity == severity]
        if not chosen:
            continue
        lines += [f"### {title}", ""]
        for category in dict.fromkeys(key[1] for key, _found in ranked if key[0] == severity):
            of_category = [finding for finding in chosen if finding.category == category]
            lines.append(f"**{category}** ({len(of_category)})")
            lines.append("")
            for finding in of_category[:LISTED]:
                where = f"`{finding.node}`: " if finding.node else ""
                lines.append(f"- {where}{finding.detail}")
            if len(of_category) > LISTED:
                lines.append(f"- and {len(of_category) - LISTED} more, in the findings table")
            lines.append("")
    return lines


def _typing_lines(typing: TypingReport) -> list[str]:
    lines = [
        "## Typing of the genomes", "",
        "Every genome of the alignment was typed with the scheme under two rules. The "
        "`single` rule is that of CanSNPer2 2.2.1: one canSNP per node, and the deepest "
        "derived node whose chain from the root holds no ancestral node, enough derived "
        "nodes, and a derived parent. The `set` rule takes the majority of "
        "each node's canSNPs and may pass a node without a decision. The reference genome "
        "and the outgroups are not counted. Percentages are of the genomes counted.", "",
        *table(SUMMARY_COLUMNS, typing.summary),
    ]

    if typing.stops:
        lines += ["### Where the single rule stopped", "",
                  "The node at which the chain of the single rule broke, or whose missing "
                  "base kept the node below it from being called, for genomes the set rule "
                  "placed deeper.", "",
                  *table(("node", "reason", "genomes"), typing.stops[:WORST])]
        if len(typing.stops) > WORST:
            lines += [f"{len(typing.stops) - WORST} further rows are in the genome table.", ""]

    if len(typing.rank_rows) > 1:
        lines += ["### The canSNPs by rank", "",
                  "Counts summed over the canSNPs of each rank and over all genomes; "
                  "`snps_disagreeing` are canSNPs that disagreed with their set or with the "
                  "path of a genome at least once.", "",
                  *table(RANK_COLUMNS, typing.rank_rows)]

    at = {column: SNP_COLUMNS.index(column) for column in SNP_COLUMNS}

    def disagreement(row: Sequence[str]) -> int:
        return max(int(row[at["discordant_in_set"]]),
                   int(row[at["ancestral_on_path"]]) + int(row[at["derived_off_path"]]))

    shown = ("node", "rank", "contig", "contig_position", "genomes", "missing",
             "discordant_in_set", "ancestral_on_path", "derived_off_path")
    disagreeing = sorted((row for row in typing.snp_rows if disagreement(row)),
                         key=lambda row: (-disagreement(row), row[at["node"]],
                                          int(row[at["rank"]])))
    if disagreeing:
        lines += ["### canSNPs that disagreed most", "",
                  f"{len(disagreeing)} of {len(typing.snp_rows)} canSNPs disagreed with the "
                  f"other canSNPs of their node or with the path of a genome at least once.", "",
                  *table(shown, [[row[at[column]] for column in shown]
                                 for row in disagreeing[:WORST]])]
    lacking = sorted((row for row in typing.snp_rows if int(row[at["missing"]])),
                     key=lambda row: (-int(row[at["missing"]]), row[at["node"]],
                                      int(row[at["rank"]])))
    if lacking:
        lines += ["### canSNPs most often without a base", "",
                  f"{len(lacking)} of {len(typing.snp_rows)} canSNPs had no base in at least "
                  f"one genome.", "",
                  *table(shown[:6], [[row[at[column]] for column in shown[:6]]
                                     for row in lacking[:WORST]])]

    if typing.unreached:
        lines += ["### Nodes these genomes do not test", "",
                  f"{len(typing.unreached)} nodes with canSNPs lie on the path of no genome: "
                  + ", ".join(f"`{node}`" for node in typing.unreached[:WORST * 2])
                  + (" and others" if len(typing.unreached) > WORST * 2 else "") + ".", ""]

    column = {name: GENOME_COLUMNS.index(name) for name in GENOME_COLUMNS}
    kept = ("genome", "single_call", "single_stop_node", "single_stop_reason", "set_call",
            "derived_off_path", "discordant")
    notable = [row for row, typed in zip(typing.genome_rows(), typing.genomes, strict=True)
               if typed.role == ISOLATE and (typed.set_call is None or typed.discordant
                                             or typed.single_call != typed.set_call)]
    if notable:
        notable.sort(key=lambda row: (row[column["set_call"]] != NA,
                                      row[column["discordant"]] != "yes",
                                      row[column["genome"]]))
        lines += ["### Genomes the two rules type differently, or not at all", "",
                  f"{len(notable)} genomes have no set call, are discordant, or get another "
                  f"call under the single rule than under the set rule.", "",
                  *table(kept, [[row[column[name]] for name in kept]
                                for row in notable[:WORST]])]
        if len(notable) > WORST:
            lines += [f"{len(notable) - WORST} further genomes are in the genome table.", ""]
    return lines


def markdown(scheme: Scheme, findings: Sequence[Finding], typing: TypingReport | None = None,
             run: Sequence[tuple[str, str]] = ()) -> str:
    """The report as Markdown text.

    Args:
        scheme: The scheme that was audited.
        findings: Its findings, in the order they are to be listed.
        typing: The typing of the genomes, when an alignment was given.
        run: What the audit was run on, as (item, value).
    """
    lines = [f"# Evaluation of the '{scheme.organism}' scheme", ""]
    if run:
        lines += ["## The run", "", *table(("item", "value"), run)]
    lines += _scheme_lines(scheme)
    lines += _finding_lines(findings)
    if typing is not None:
        lines += _typing_lines(typing)
    return "\n".join(lines).rstrip("\n") + "\n"


def write_markdown(path: str, scheme: Scheme, findings: Sequence[Finding],
                   typing: TypingReport | None = None,
                   run: Sequence[tuple[str, str]] = ()) -> None:
    with open(path, "w", newline="\n") as handle:
        handle.write(markdown(scheme, findings, typing, run))
