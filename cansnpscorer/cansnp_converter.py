import argparse
from cansnpscorer.helpers.file_utils import convert_data, save_to_file, recognize_hc_file_format

def main(infile: str, prefix: str, format_in: str, format_out: str):
    input_format = format_in or recognize_file_format(infile)
    output_format = format_out

    with open(infile, 'r') as f:
        data = f.read()

    converted_data = convert_data(data, input_format, output_format)

    print(converted_data)
    save_to_file(converted_data, prefix, output_format)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert hierarchical data formats')
    parser.add_argument("-i","--infile", type=str, help="The input file name")
    parser.add_argument("-p","--prefix", type=str, help="The output file name")
    parser.add_argument("--format_in", type=str, default=None, help="The format of the input file ('yaml', 'json', 'newick', 'child_parent')")
    parser.add_argument("--format_out", type=str, required=True, help="The format of the output file ('yaml', 'json', 'newick', 'child_parent')")

    args = parser.parse_args()
    main(args.infile, args.prefix, args.format_in, args.format_out)

