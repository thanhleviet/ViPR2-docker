#!/usr/bin/env python

# Filtering overlapped POS in VCF file. Criteria: Select POS with highest QUAL values
# Thanh Le Viet
# Oct 2017 thanhlv@oucru.org

import vcf
import argparse
from collections import  Counter
import logging
import sys

LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

def cmdline_parser():
    """
    creates an argparse instance
    """

    parser = argparse.ArgumentParser(__doc__)

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Optional: be verbose")
    parser.add_argument("--force",
                        action="store_true",
                        dest="force_overwrite",
                        help="Optional: force overwriting of files")
    parser.add_argument("--debug",
                        action="store_true",
                        dest="debug",
                        help="Optional: enable debugging")
    parser.add_argument("-i", "--vcf-in",
                        dest="vcf_in",
                        help="VCF input file")
    parser.add_argument("-o", "--vcf-out",
                        dest="vcf_out",
                        help="VCF output file")
    return parser

def dup_index(List, value):
    """
    Find index of duplicated elements in a list
    """
    return [i for i,x in enumerate(List) if x == value]

def write_vcf(out_vcf, template, vcf_record):
    """
    Write to a vcf file filtered VCF record
    """
    out_vcf_ = open(out_vcf, 'w')
    vcf_writer = vcf.Writer(out_vcf_, template)
    for f in vcf_record:
        vcf_writer.write_record(f)
    out_vcf_.close()

def vcf_filter(in_vcf, out_vcf):
    vcf_recs = [x.POS for x in vcf.Reader(open(in_vcf, "r"))]
    vcf_ = [x for x in vcf.Reader(open(in_vcf, "r"))]

    #Extract overlapped POS value
    overlap = [k for k, v in Counter(vcf_recs).items() if v > 1]
    if len(overlap) > 0:
        LOG.info("Found {} overlapped!".format(len(overlap)))
        overlap.sort()
        #Find index for each overlapped POS value
        dup_pos = dict((x, dup_index(vcf_recs, x)) for x in set(vcf_recs) if vcf_recs.count(x) > 1)

        drop = []
        for x in overlap:
            QUAL = [vcf_[i].QUAL for i in dup_pos[x]]
            MAX_QUAL = max(QUAL)
            for i, k in enumerate(QUAL):
                if k != MAX_QUAL:
                    drop.append(dup_pos[x][i])

        filtered = [x for i,x in enumerate(vcf_) if i not in drop]
        dropped = [x for i,x in enumerate(vcf_) if i in drop]

        template = vcf.Reader(open(in_vcf, "r"))

        write_vcf(out_file, template, filtered)

        dropped_vcf = "exc_dup_{}".format(out_file)
        write_vcf(dropped_vcf, template, dropped)
    else:
        LOG.info("No overlapped POS")

if __name__ == "__main__":
    parser = cmdline_parser()
    args = parser.parse_args()

    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)

    for (filename, descr, direction, mandatory) in [
            (args.vcf_in, 'VCF input file', 'in', True),
            (args.vcf_out, 'VCF output file', 'out', True)]:

        if not mandatory and not filename:
            continue

        if not filename:
            parser.error("%s argument missing." % descr)
            sys.exit(1)

    in_file = args.vcf_in
    out_file = args.vcf_out
    vcf_filter(in_vcf=in_file, out_vcf=out_file)
    LOG.info("Successful exit")