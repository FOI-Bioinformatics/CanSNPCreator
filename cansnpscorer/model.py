"""The reference genome, and the map from alignment columns onto it.

The tools work with two coordinate systems. A whole genome alignment is
addressed by column, counted over every contig of the reference and including
the gaps the aligner introduced. The reference genome is addressed by contig
and position within that contig, which is what a published canSNP scheme
records and what CanSNPer2 expects.

``Reference`` holds the contigs and converts between a position in the
concatenated genome and a contig with a local position. ``AlignmentMap``
converts an alignment column into a position in the reference genome, and
refuses to do so unless the alignment and the reference agree base for base.
"""

from array import array
from bisect import bisect_right
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

GAP_CHARACTERS = frozenset("-.~")
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


class ReferenceError(ValueError):
    """The reference genome, or something addressed against it, is unusable."""


def complement(base: str) -> str:
    """The complement of a base, upper case; anything unknown is returned as N."""
    return COMPLEMENT.get(base.upper(), "N")


def reverse_complement(sequence: str) -> str:
    return "".join(complement(base) for base in reversed(sequence))


def ungap(sequence: str) -> str:
    return "".join(base for base in sequence if base not in GAP_CHARACTERS)


@dataclass(frozen=True)
class Contig:
    """One sequence of the reference, with its offset in the concatenated genome."""

    name: str
    sequence: str
    start: int

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def end(self) -> int:
        """One past the last position of this contig, in global coordinates."""
        return self.start + self.length

    def base(self, local_position: int) -> str:
        """The base at a 0-based position within this contig."""
        if not 0 <= local_position < self.length:
            raise ReferenceError(
                f"position {local_position + 1} is outside contig '{self.name}' "
                f"of {self.length} bases"
            )
        return self.sequence[local_position]

    def window(self, local_position: int, distance: int) -> str:
        """The sequence around a position, clamped to the ends of the contig."""
        start = max(local_position - distance, 0)
        return self.sequence[start:local_position + distance]


class Reference:
    """The contigs of a reference genome, in file order."""

    def __init__(self, contigs: Sequence[Contig]) -> None:
        if not contigs:
            raise ReferenceError("the reference holds no sequence")
        self.contigs: List[Contig] = list(contigs)
        self._starts = [contig.start for contig in self.contigs]
        self._by_name: Dict[str, Contig] = {contig.name: contig for contig in self.contigs}

    @classmethod
    def from_records(cls, records: Iterable[Tuple[str, str]]) -> "Reference":
        """Build a reference from (name, sequence) pairs in file order."""
        contigs: List[Contig] = []
        start = 0
        for name, sequence in records:
            if name in {contig.name for contig in contigs}:
                raise ReferenceError(f"the reference names contig '{name}' more than once")
            if not sequence:
                raise ReferenceError(f"contig '{name}' holds no sequence")
            contigs.append(Contig(name=name, sequence=sequence, start=start))
            start += len(sequence)
        return cls(contigs)

    @property
    def length(self) -> int:
        """The total number of bases over all contigs."""
        return self.contigs[-1].end

    @property
    def sequence(self) -> str:
        """Every contig concatenated, in file order."""
        return "".join(contig.sequence for contig in self.contigs)

    @property
    def names(self) -> List[str]:
        return [contig.name for contig in self.contigs]

    def contig(self, name: str) -> Contig:
        if name not in self._by_name:
            raise ReferenceError(
                f"the reference has no contig '{name}'; it has {', '.join(self.names)}"
            )
        return self._by_name[name]

    def locate(self, position: int) -> Tuple[Contig, int]:
        """Split a 0-based global position into a contig and a 0-based local position."""
        if not 0 <= position < self.length:
            raise ReferenceError(
                f"position {position + 1} is outside the reference of {self.length} bases"
            )
        contig = self.contigs[bisect_right(self._starts, position) - 1]
        return contig, position - contig.start

    def base(self, position: int) -> str:
        contig, local_position = self.locate(position)
        return contig.base(local_position)

    def __len__(self) -> int:
        return len(self.contigs)

    def __iter__(self):
        return iter(self.contigs)


class AlignmentMap:
    """Alignment columns mapped onto positions in the reference genome.

    Two layouts of whole genome alignment are accepted, and told apart from the
    reference row itself:

    - **gapped**: the reference row carries the gaps the aligner introduced, so
      a column is not a reference position, and the row with its gaps removed
      is the whole reference genome. progressiveMauve's own output, turned into
      FASTA column by column, has this layout.
    - **projected**: every row has been projected onto reference coordinates,
      so column i is position i of the concatenated reference. Columns where
      the reference was not aligned at all are gaps in the reference row too,
      and the row may stop short of the reference's end. The x2fa conversion
      used by CanSNPer, RepGenR and CRSP, and the whole genome alignments of
      snippy and RepGenR2's simple SNP typer, have this layout.

    Either way, building the map checks the alignment against the reference
    base for base, since a mismatch would move every position reported.
    """

    def __init__(self, aligned_reference: str, reference: Reference, examples: int = 3) -> None:
        self.reference = reference
        expected = reference.sequence.upper()
        observed = ungap(aligned_reference).upper()

        if len(observed) == len(expected):
            self.projected = False
            positions = array("q")
            position = 0
            for base in aligned_reference:
                if base in GAP_CHARACTERS:
                    positions.append(-1)
                else:
                    positions.append(position)
                    position += 1
            mismatches = [
                index for index, (left, right) in enumerate(zip(observed, expected))
                if left != right
            ]
        elif len(aligned_reference) <= len(expected):
            self.projected = True
            positions = array("q", range(len(aligned_reference)))
            mismatches = []
            for index, base in enumerate(aligned_reference):
                if base in GAP_CHARACTERS:
                    positions[index] = -1
                elif base.upper() != expected[index]:
                    mismatches.append(index)
            observed = aligned_reference.upper()
        else:
            raise ReferenceError(
                f"the alignment's reference row holds {len(observed)} bases once its gaps "
                f"are removed, but the reference genome holds {len(expected)}. The first "
                f"record of the alignment must be the reference genome, and the reference file "
                f"must hold the contigs the alignment concatenates, in that order"
            )

        if mismatches:
            shown = []
            for index in mismatches[:examples]:
                contig, local_position = reference.locate(index)
                shown.append(
                    f"{contig.name}:{local_position + 1} has {observed[index]} in the "
                    f"alignment and {expected[index]} in the reference"
                )
            raise ReferenceError(
                f"the alignment's reference row differs from the reference genome at "
                f"{len(mismatches)} positions: " + "; ".join(shown)
            )

        self._positions = positions
        self._columns: Optional[array] = None

    @property
    def columns(self) -> int:
        return len(self._positions)

    @property
    def unaligned(self) -> int:
        """Reference positions that no column of the alignment covers."""
        covered = sum(1 for position in self._positions if position >= 0)
        return self.reference.length - covered

    def to_reference(self, column: int) -> Optional[int]:
        """The 0-based global position for a column, or None if the reference is gapped there."""
        if not 0 <= column < self.columns:
            raise ReferenceError(f"column {column} is outside the alignment")
        position = self._positions[column]
        return None if position < 0 else position

    def to_column(self, position: int) -> Optional[int]:
        """The column holding a 0-based global position, or None if no column does."""
        if not 0 <= position < self.reference.length:
            raise ReferenceError(
                f"position {position + 1} is outside the reference of "
                f"{self.reference.length} bases")
        if self._columns is None:
            columns = array("q", [-1]) * self.reference.length
            for column, found in enumerate(self._positions):
                if found >= 0:
                    columns[found] = column
            self._columns = columns
        column = self._columns[position]
        return None if column < 0 else column

    def locate(self, column: int) -> Optional[Tuple[Contig, int]]:
        """The contig and 1-based position for a column, or None if the reference is gapped."""
        position = self.to_reference(column)
        if position is None:
            return None
        contig, local_position = self.reference.locate(position)
        return contig, local_position + 1
