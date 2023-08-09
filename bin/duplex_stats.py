#!/usr/bin/env python
"""Count duplex and simplex reads in xam file."""
import argparse

import pysam


def main():
    """Run entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument('xam')
    parser.add_argument('outname')
    args = parser.parse_args()

    # Prepare input files
    xam = pysam.AlignmentFile(args.xam, check_sq=False)

    # Count simplex and duplex reads in a xam
    sx = dx = 0
    for read in xam.fetch(until_eof=True):
        if read.get_tag('dx') == 1:
            dx += 1
        else:
            sx += 1

    # Save counts to output
    with open(f'{args.outname}', 'w') as out:
        out.write("Filename,Duplex,Paired,Simplex\n")
        out.write(f"{args.xam},{dx},{dx*2},{sx}\n")


if __name__ == '__main__':
    main()
