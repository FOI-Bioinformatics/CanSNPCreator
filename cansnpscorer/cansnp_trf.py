"""Running Tandem Repeats Finder over a reference genome, and writing its
repeats as BED."""

import logging
import os
import shutil
import subprocess
from typing import List, Tuple

TRF_EXECUTABLE = "trf"

#: (sequence ID, 1-based start, 1-based end, period, copy number) per repeat.
Repeat = Tuple[str, int, int, int, str]


def parse_ngs(text: str) -> List[Repeat]:
    """
    Read the repeats of Tandem Repeats Finder's -ngs output.

    Each sequence begins with a line '@' followed by its header, whose first
    word is the sequence ID, as FASTA defines it. Each repeat line then gives
    its start, end, period and copy number first, 1-based and inclusive.

    Raises:
        ValueError: If a repeat line comes before any sequence, or cannot be read.
    """
    repeats: List[Repeat] = []
    sequence = None
    for number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        if line.startswith("@"):
            words = line[1:].split()
            sequence = words[0] if words else None
            continue
        if sequence is None:
            raise ValueError(f"line {number} of the TRF output comes before any sequence")
        fields = line.split()
        try:
            start, end, period = int(fields[0]), int(fields[1]), int(fields[2])
            copies = fields[3]
        except (IndexError, ValueError):
            raise ValueError(f"line {number} of the TRF output is not a repeat: {line}") from None
        repeats.append((sequence, min(start, end), max(start, end), period, copies))
    return repeats


def write_bed(repeats: List[Repeat], path: str) -> None:
    """Write repeats as BED: 0-based, half open, in the order TRF met the sequences."""
    order = {}
    for sequence, *_rest in repeats:
        order.setdefault(sequence, len(order))
    with open(path, "w", newline="\n") as handle:
        for sequence, start, end, period, copies in sorted(
                repeats, key=lambda repeat: (order[repeat[0]], repeat[1], repeat[2])):
            handle.write(f"{sequence}\t{start - 1}\t{end}\tperiod={period},copies={copies}\n")


def run_trf(reference: str, output: str, match=2, mismatch=7, delta=7, pm=80, pi=10,
            minscore=10, max_period=2000) -> int:
    """
    Run Tandem Repeats Finder over a reference and write its repeats as BED.

    Args:
        reference: The reference genome in FASTA format.
        output: Where to write the BED file.
        match, mismatch, delta, pm, pi, minscore, max_period: TRF's parameters,
            in the order TRF takes them.

    Returns:
        A process exit status: 0 when the BED file was written, 1 otherwise.
        Nothing is written when TRF fails, since an empty file would read as a
        genome without repeats.
    """
    if not os.path.isfile(reference):
        logging.error(f"the reference {reference} does not exist")
        return 1
    if shutil.which(TRF_EXECUTABLE) is None:
        logging.error(
            f"'{TRF_EXECUTABLE}' is not on the PATH. Install Tandem Repeats Finder, for "
            "example with 'conda install -c bioconda trf'.")
        return 1

    command = [
        TRF_EXECUTABLE, reference,
        str(match), str(mismatch), str(delta),
        str(pm), str(pi), str(minscore), str(max_period),
        "-h", "-ngs",
    ]
    logging.info(f"Running {' '.join(command)}")

    try:
        completed = subprocess.run(command, capture_output=True, text=True)
    except OSError as error:
        logging.error(f"Could not run Tandem Repeats Finder: {error}")
        return 1

    if completed.returncode != 0:
        logging.error(
            f"Tandem Repeats Finder exited with status {completed.returncode}, so nothing was "
            f"written. {completed.stderr.strip()}")
        return 1

    try:
        repeats = parse_ngs(completed.stdout)
    except ValueError as error:
        logging.error(error)
        return 1
    if not repeats:
        logging.warning(
            "Tandem Repeats Finder reported no repeat at all, which is unusual for a bacterial "
            "genome; check its parameters before relying on the empty file")
    write_bed(repeats, output)
    logging.info(f"Wrote {len(repeats)} tandem repeats to {output}")
    return 0
