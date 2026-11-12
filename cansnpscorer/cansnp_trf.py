"""Running Tandem Repeats Finder over a reference genome."""

import logging
import shutil
import subprocess

TRF_EXECUTABLE = "trf"


def run_trf(input_file: str, output_file: str, match=2, mismatch=7, delta=7, pm=80, pi=10,
            minscore=10, max_period=2000) -> int:
    """
    Run Tandem Repeats Finder and keep its report.

    Args:
        input_file: Path to the reference genome in FASTA format.
        output_file: Where to write the repeat report.
        match: Matching weight. Default is 2.
        mismatch: Mismatching penalty. Default is 7.
        delta: Indel penalty. Default is 7.
        pm: Match probability as a whole number. Default is 80.
        pi: Indel probability as a whole number. Default is 10.
        minscore: Minimum alignment score to report. Default is 10.
        max_period: Maximum period size to report. Default is 2000.

    Returns:
        A process exit status: 0 when the report was written, 1 otherwise.
        A failed run used to leave an empty report behind, which the evaluator
        then read as a genome with no repeats anywhere.
    """
    if shutil.which(TRF_EXECUTABLE) is None:
        logging.error(
            f"'{TRF_EXECUTABLE}' is not on the PATH. Install Tandem Repeats Finder, for "
            "example with 'conda install -c bioconda trf'.")
        return 1

    command = [
        TRF_EXECUTABLE, input_file,
        str(match), str(mismatch), str(delta),
        str(pm), str(pi), str(minscore), str(max_period),
        "-h", "-ngs",
    ]
    logging.info(f"Running {' '.join(command)}")

    try:
        with open(output_file, "w") as handle:
            completed = subprocess.run(command, stdout=handle, stderr=subprocess.PIPE, text=True)
    except OSError as error:
        logging.error(f"Could not run Tandem Repeats Finder: {error}")
        return 1

    if completed.returncode != 0:
        logging.error(
            f"Tandem Repeats Finder exited with status {completed.returncode}; "
            f"{output_file} is incomplete and must not be used. "
            f"{completed.stderr.strip()}")
        return 1

    logging.info(f"Wrote the tandem repeat report to {output_file}")
    return 0
