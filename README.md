# CanSNPCreator

## Introduction

CanSNPCreator, formerly known as CanSNPScorer, is an advanced tool designed for the extraction, evaluation, and analysis of Canonical Single Nucleotide Polymorphisms (canSNPs) in microbial genomes. With enhanced support for various hierarchical clustering file formats, CanSNPCreator provides a robust platform for genomic analysis.

## Features

- **SNP Extraction**: Extract SNPs from genomic sequence alignments.
- **canSNP Extraction**: Identify canSNPs from multiple hierarchical clustering file formats including YAML and repgenr HC formats.
- **Evaluation**: Assess canSNPs with additional genomic context using a comprehensive set of parameters.
- **TRF Integration**: Incorporate Tandem Repeats Finder results for more nuanced canSNP evaluation.


## Requirements

- Python 3.9 or higher
- Additional Python libraries as specified in `requirements.txt`
- Tandem Repeat Finder (TRF)

## Installation

1. Clone the repository:

    ```bash
    git clone https://github.com/yourusername/CanSNPCreator.git
    ```

2. Navigate to the cloned directory:

    ```bash
    cd CanSNPCreator
    ```

3. Install the required Python packages:

    ```bash
    pip install -r requirements.txt
    ```

## Usage

### Basic Commands

1. **Extract canSNPs**

    ```bash
    CanSNPCreator extract --csvfile <path-to-csvfile> --hcfile <path-to-hcfile> --fastafile <path-to-fastafile> -o <output-file>
    ```

3. **Perform TRF on reference genome**

    ```bash
    CanSNPCreator trf --input <input-file> --output <output-file>
    ```


2. **Evaluate canSNPs**

    ```bash
    CanSNPCreator evaluate --snp-file <path-to-snpfile> --reference-file <path-to-reference-file> --gff-file <path-to-gff-file> --trf-file <path-to-trf-file> --output-file <output-file>
    ```




### Advanced Options

For a complete list of options, run:

```bash
CanSNPCreator -h
```

## Contributing

If you'd like to contribute, please fork the repository and create a pull request, or simply open an issue for any bug reports or feature requests.

## License

This project is licensed under the PGL3 License. See the `LICENSE` file for details.

## Contact

For any questions or suggestions, feel free to [open an issue](https://github.com/Bioinformatics/CanSNPScorer/issues) or contact the main contributors.
