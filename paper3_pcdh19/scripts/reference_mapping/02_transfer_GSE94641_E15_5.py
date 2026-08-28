#!/usr/bin/env python3
"""Run E15.5-focused GSE94641 reference-PCA/kNN annotation transfer."""

from __future__ import print_function

import argparse
import json
import os
import sys

from gse94641_label_transfer import GSE94641LabelTransferWorkflow, MappingError


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", required=True)
    parser.add_argument("--reference-root", required=True)
    parser.add_argument("--validation-root", required=True)
    parser.add_argument("--query-root", required=True)
    parser.add_argument("--sample-key", required=True)
    parser.add_argument("--output-root", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    with open(args.lock, "r") as handle:
        configuration = json.load(handle)
    paths = {
        "reference_root": os.path.abspath(args.reference_root),
        "validation_root": os.path.abspath(args.validation_root),
        "query_root": os.path.abspath(args.query_root),
        "sample_key": os.path.abspath(args.sample_key),
        "output_root": os.path.abspath(args.output_root),
    }
    GSE94641LabelTransferWorkflow(configuration, paths).run()


if __name__ == "__main__":
    try:
        main()
    except MappingError as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        sys.exit(2)
