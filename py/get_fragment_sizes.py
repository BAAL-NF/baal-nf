"""
Get the fragment sizes extracted by phantompeakqualtools
and pick the highest-probability nonzero result. 

Exits with status 3 if no suitable numbers are found
"""

from argparse import ArgumentParser, FileType
from csv import reader


def get_input():
    parser = ArgumentParser()
    parser.add_argument(
        "input", help="input file produced by run_spp.R", type=FileType("r")
    )
    parser.add_argument(
        "output", help="name of file to output result to", type=FileType("w")
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = get_input()
    spp = reader(args.input, delimiter="\t")

    # The fragment lengths are given as a comma-separated string
    fragment_lengths = map(lambda x: int(x), next(spp)[2].split(","))

    for fragment_length in fragment_lengths:
        if fragment_length != 0:
            # write the fragment length to a file, ignoring the sign
            args.output.write(f"{abs(fragment_length)}\n")
            exit(0)

    # no nonzero fragment lengths provided. Exit with an error.
    exit(3)
