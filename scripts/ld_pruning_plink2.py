#!/usr/bin/env python
import sys
import argparse
import subprocess


def main(argv):
    parser = argparse.ArgumentParser(description='Perform LD pruning using plink2')
    parser.add_argument('--vcf', help='Input vcf. Multiple files can be specified separated by a space.', nargs='+',
                        required=True)
    parser.add_argument('--output_base', help='Output bases. Multiple files can be specified separated by a space.',
                        nargs='+', required=True)
    parser.add_argument('--sex', help='Individual sex mapping files. Multiple files can be specified separated '
                                      'by a space.', nargs='+', required=True)
    parser.add_argument('--threads', help='Number of CPUs [1]', default=1)
    args = parser.parse_args()
    for i, o, s in zip(args.vcf, args.output_base, args.sex):
        subprocess.call(f"plink2 --threads {args.threads} --allow-extra-chr --vcf {i} --indep-pairwise 50 kb 1 0.1 "
                        f"--out {o} --update-sex {s}", shell=True)


if __name__ == '__main__':
    main(sys.argv[1:])
