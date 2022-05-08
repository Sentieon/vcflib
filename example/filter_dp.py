#!/usr/bin/env python

"""
An example script that utilizes vcflib to remove variants in a VCF file by DP
"""

from __future__ import print_function

import argparse
import os
import os.path
import sys

import vcflib


def filter_vcf(in_vcf, out_vcf, min_dp=10):
    """ Filter variants in the input VCF with a INFO/DP < min_dp """
    for variant in in_vcf:
        dp = variant.info.get("DP", None)  # INFO fields are stored as a dict
        if dp is not None and dp >= min_dp:
            out_vcf.emit(variant)
    return


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input_vcf", required=True, help="The input VCF")
    parser.add_argument("--output_vcf", required=True, help="The output VCF")
    parser.add_argument(
        "--min_dp", default=10, type=int, help="The minimum DP for a variant"
    )
    parser.add_argument(
        "--n_threads",
        default=1,
        type=int,
        help="Number of threads to use for VCF processing",
    )
    parser.add_argument(
        "--step_size",
        default=10 * 1000 * 1000,
        type=int,
        help="Size of the shards",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    assert os.path.isfile(args.input_vcf)
    in_vcf = vcflib.VCF(str(args.input_vcf))
    # `VCF()` accepts an option mode argument. Behavior is similar to mode in `open()`
    out_vcf = vcflib.VCF(str(args.output_vcf), "wb")
    out_vcf.copy_header(
        in_vcf
    )  # Copy headers from the input VCF to the output VCF
    out_vcf.emit_header()  # Write the VCF header in the output VCF

    if args.n_threads < 2:
        # Run the filtering in the current process without parallelization
        filter_vcf(in_vcf, out_vcf, args.min_dp)
    else:
        ## Run the filtering in multiple processes - requires a VCF with contig
        ## lengths in the header
        # Initialize the Sharder object with the number of threads
        sharder = vcflib.Sharder(args.n_threads)
        # Partition the genome into equal-sized shards
        try:
            contig_lengths = [
                (contig, 0, int(d["length"]))
                for contig, d in in_vcf.contigs.items()
            ]
        except KeyError:
            return (
                "This example script requires a VCF with contig lengths in"
                " the header when using multiple threads"
            )
        shards = sharder.cut(contig_lengths, args.step_size)
        # Run all of the shards in parallel
        sharder.run(shards, filter_vcf, [], in_vcf, out_vcf, args.min_dp)

    out_vcf.close()  # Close the output VCF and write the index file
    return os.EX_OK


if __name__ == "__main__":
    sys.exit(main())
