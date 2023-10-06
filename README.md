# CanSNPScorer

## Introduction

CanSNPScorer is a comprehensive tool designed for the extraction and evaluation of Canonical Single Nucleotide Polymorphisms (canSNPs) in microbial genomes. The tool integrates multiple genomic data types and formats, offering a unified platform for canSNP analysis.

## Features

- **SNP Extraction**: Extract SNPs from genomic sequence alignments.
- **canSNP Extraction**: Identify canSNPs based on Hierarchical Clustering (HC).
- **Evaluation**: Evaluate canSNPs with additional genomic context.
- **TRF**: Tandem Repeats Finder integration for canSNP evaluation.

## Requirements

- Python 3.9 or higher
- Additional Python libraries as specified in `requirements.txt`
- Tandem Repeat Finder (trf)

## Installation

1. Clone the repository:

    ```bash
    git clone https://github.com/yourusername/CanSNPScorer.git
    ```

2. Navigate to the cloned directory:

    ```bash
    cd CanSNPScorer
    ```

3. Install the required Python packages:

    ```bash
    pip install -r requirements.txt
    ```

## Usage

### Basic Commands

1. **Extract canSNPs**

    ```bash
    CanSNPScorer extract --csvfile <path-to-csvfile> --hcfile <path-to-hcfile> --fastafile <path-to-fastafile> -o <output-file>
    ```

2. **Evaluate canSNPs**

    ```bash
    CanSNPScorer evaluate --snp-file <path-to-snpfile> --reference-file <path-to-reference-file> --gff-file <path-to-gff-file> --trf-file <path-to-trf-file> --output-file <output-file>
    ```

3. **Perform TRF on canSNPs**

    ```bash
    CanSNPScorer trf --input <input-file> --output <output-file>
    ```

### Advanced Options

For a complete list of options, run:

```bash
CanSNPScorer -h
```

## Contributing

If you'd like to contribute, please fork the repository and create a pull request, or simply open an issue for any bug reports or feature requests.

## License

This project is licensed under the PGL3 License. See the `LICENSE` file for details.

## Contact

For any questions or suggestions, feel free to [open an issue](https://github.com/Bioinformatics/CanSNPScorer/issues) or contact the main contributors.
