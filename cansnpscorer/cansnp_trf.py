import subprocess
import logging

def run_trf(input_file: str, output_file: str, match=2, mismatch=7, delta=7, pm=80, pi=10, minscore=10, max_period=2000):
    """
    Run the Tandem Repeats Finder (TRF) tool with specific parameters.

    Parameters:
    - input_file (str): Path to the input FASTA file.
    - output_file (str): Path to the output file where TRF results will be stored.
    - match (int): Matching weight. Default is 2.
    - mismatch (int): Mismatching penalty. Default is 7.
    - delta (int): Indel penalty. Default is 7.
    - pm (int): Match probability (whole number). Default is 80.
    - pi (int): Indel probability (whole number). Default is 10.
    - minscore (int): Minimum alignment score to report. Default is 10.
    - max_period (int): Maximum period size to report. Default is 2000.

    Returns:
    - None
    """
    try:
        logging.info(f"Running TRF on {input_file}...")

        # Define the TRF command with parameters
        trf_command = [
            "trf",
            input_file,
            str(match),
            str(mismatch),
            str(delta),
            str(pm),
            str(pi),
            str(minscore),
            str(max_period),
            "-h",
            "-ngs"
        ]

        # Run the command and save the output
        with open(output_file, "w") as f_out:
            subprocess.run(trf_command, stdout=f_out)

        logging.info(f"TRF analysis completed. Results are saved in {output_file}.")

    except Exception as e:
        logging.error(f"An error occurred while running TRF: {e}")

def main():
    parser = argparse.ArgumentParser(description='Run the Tandem Repeats Finder (TRF) tool with specific parameters.')
    parser.add_argument('--input', required=True, help='Input fasta file for TRF')
    parser.add_argument('--output', required=True, help='Output file for TRF results')
    parser.add_argument('--match', type=int, default=2, help='Matching weight. Default is 2.')
    parser.add_argument('--mismatch', type=int, default=7, help='Mismatching penalty. Default is 7.')
    parser.add_argument('--delta', type=int, default=7, help='Indel penalty. Default is 7.')
    parser.add_argument('--pm', type=int, default=80, help='Match probability (whole number). Default is 80.')
    parser.add_argument('--pi', type=int, default=10, help='Indel probability (whole number). Default is 10.')
    parser.add_argument('--minscore', type=int, default=10, help='Minimum alignment score to report. Default is 10.')
    parser.add_argument('--max_period', type=int, default=2000, help='Maximum period size to report. Default is 2000.')

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info("Starting TRF analysis...")


    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info("Starting TRF analysis...")

    run_trf(
        trf_input,
        trf_output,
        trf_match,
        trf_mismatch,
        trf_delta,
        trf_pm,
        trf_pi,
        trf_minscore,
        trf_max_period
    )

if __name__ == '__main__':
    main()
