#!/usr/bin/env python
# -*- coding: utf-8 -*-

# From: https://github.com/brentp/bio-playground
# Refer to included "bio-playground_LICENSE.txt" file for
#     license and copyright statement.

"""
Guess the (quality-score) encoding of a FASTQ file.
"""

from __future__ import with_statement, division, print_function

import argparse
import fileinput
import itertools
import operator
import sys

from collections import Counter

#  Note that the theoretical maximum for all encodings is 126.
#  The upper limits below are for "typical" data only.
RANGES = {
    'Sanger': (33, 73),
    'Illumina-1.8': (33, 74),
    'Solexa': (59, 104),
    'Illumina-1.3': (64, 104),
    'Illumina-1.5': (66, 105)
}

# The threshold to decide between Illumina-1.3 and Illumina-1.5
# based upon how common "B" is. The threshold insists it is
# within the Nth most common quality scores.
# N.B. needs to be conservative, as this is applied per input line.
N_MOST_COMMON_THRESH = 4


def get_qual_range(qual_str):
    """
    >>> get_qual_range("DLXYXXRXWYYTPMLUUQWTXTRSXSWMDMTRNDNSMJFJFFRMV")
    (68, 89...)
    """

    qual_val_counts = Counter([ord(qual_char) for qual_char in qual_str])

    min_base_qual = min(qual_val_counts.elements())
    max_base_qual = max(qual_val_counts.elements())

    return (min_base_qual, max_base_qual, qual_val_counts)


def get_encodings_in_range(rmin, rmax, ranges=RANGES):
    valid_encodings = []
    for encoding, (emin, emax) in ranges.items():
        if rmin >= emin and rmax <= emax:
            valid_encodings.append(encoding)
    return valid_encodings


def heuristic_filter(valid, qual_val_counts,
                     disable_uncertain_heuristics=False):
    """Apply heuristics to particular ASCII value scores
       to try to narrow-down the encoding, beyond min/max.
    """

    is_valid_modified = False

    if len(valid) > 1 and 'Illumina-1.5' in valid:
        # 64–65: Phread+64 quality scores of 0–1 ('@'–'A')
        #        unused in Illumina 1.5+
        if qual_val_counts[64] > 0 or qual_val_counts[65] > 0:
            valid.remove('Illumina-1.5')

            return valid, True

        if not disable_uncertain_heuristics:
            # 66: Phread+64 quality score of 2 'B'
            #     used by Illumina 1.5+ as QC indicator
            if 66 in map(operator.itemgetter(0),
                         qual_val_counts.most_common(N_MOST_COMMON_THRESH)):
                print("# A large number of 'B' scores (value 2, ASCII 66) "
                      "were detected, making it likely that this encoding is "
                      "Illumina-1.5, which is now the only guess made.",
                      file=sys.stderr)

                valid = ['Illumina-1.5']

                return valid, True

    return valid, is_valid_modified


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.
                                     RawDescriptionHelpFormatter)
    parser.add_argument('input_FASTQ', nargs='?', default='-',
                        help="The full path to the "
                        "FASTQ to process. or '-' or nothing to process "
                        "from STDIN. Input files can be gzipped.")
    parser.add_argument("-n", "--num_lines_to_test",
                        help="Number of qual lines to test. "
                        "By default, tests until the end of file or until "
                        "it it possible to determine a unique encoding.",
                        type=int, default=-1)
    parser.add_argument('-d', '--disable_early_stop_from_heuristics',
                        action='store_true', help="Disable any early-stopping "
                        "from heuristics beyond min/max. Setting this "
                        "prevents lack of information on min/max, "
                        "due to not parsing most of the file.")
    parser.add_argument('-D', '--disable_uncertain_heuristics',
                        action='store_true', help="Disable any heuristics "
                        "that involve a non-zero probability of error. "
                        "Does not disable the use of typical ranges, "
                        "which are themselves not entirely reliable.")

    args = parser.parse_args()

    gmin = 99
    gmax = 0
    valid = []
    stop_updating_valid = False

    err_exit = False

    input_file = fileinput.input(args.input_FASTQ,
                                 openhook=fileinput.hook_compressed)

    # every 4th line following the 4th line
    quality_scores = itertools.islice(input_file, 3, None, 4)

    for i, line in enumerate(quality_scores):
        if i == 0:
            input_filename_for_disp = fileinput.filename()

            if fileinput.isstdin():
                input_filename_for_disp = 'STDIN'

            print("# reading qualities from "
                  "{}".format(input_filename_for_disp), file=sys.stderr)

        lmin, lmax, qual_val_counts = get_qual_range(line.rstrip())

        if lmin < gmin or lmax > gmax:
            gmin, gmax = min(lmin, gmin), max(lmax, gmax)

            if not stop_updating_valid:
                valid = get_encodings_in_range(gmin, gmax)

                valid, is_valid_modified = \
                    heuristic_filter(valid, qual_val_counts,
                                     args.disable_uncertain_heuristics)

            if len(valid) == 0:
                print("no encodings for range: "
                      "{}".format((gmin, gmax)), file=sys.stderr)
                err_exit = True
                break

            if len(valid) == 1 and args.num_lines_to_test == -1:
                # parsed entire file and found unique guess
                if (not
                    (is_valid_modified and args.
                     disable_early_stop_from_heuristics)):
                    # do not break if the unique guess is from a heuristic
                    # and early stopping due to heuristics is prohibited
                    break
                elif is_valid_modified:
                    stop_updating_valid = True

        if args.num_lines_to_test > 0 and i > args.num_lines_to_test:
            # parsed up to specified portion; return current guess(es)
            break

    input_file.close()

    if err_exit:
        sys.exit(1)
    else:
        print("{}\t{}\t{}".format(valid, gmin, gmax),
              file=sys.stderr)


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |
                       doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
