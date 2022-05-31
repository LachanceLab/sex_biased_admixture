#!/usr/bin/env python
import sys
import argparse
import subprocess


def main(argv):
    parser = argparse.ArgumentParser(description='Convert vcf to plink file')
    parser.add_argument('--vcf', help='Input vcf. Multiple files can be specified separated by a space.', nargs='+',
                        required=True)
    parser.add_argument('--output_base', help='Output bases. Multiple files can be specified separated by a space.',
                        nargs='+', required=True)
    parser.add_argument('--sex', help='Individual sex mapping files. Multiple files can be specified separated '
                                      'by a space.', nargs='+', required=True)
    parser.add_argument('--pruned', help='Plink *prune.in files. Multiple files can be specified separated by a space.',
                        nargs='+', default=[])
    parser.add_argument('--threads', help='Number of CPUs [1]', default=1)
    args = parser.parse_args()
    if len(args.pruned) > 0:
        for i, o, s, p in zip(args.vcf, args.output_base, args.sex, args.pruned):
            subprocess.call(f"plink2 --threads {args.threads} --chr-set 40 --allow-extra-chr --vcf {i} --make-bed "
                            f"--extract {p} --out {o} --pca 20 approx --update-sex {s} --max-alleles 2 --maf 0.01",
                            shell=True)
    else:
        for i, o, s in zip(args.vcf, args.output_base, args.sex):
            subprocess.call(f"plink2 --threads {args.threads} --chr-set 40 --allow-extra-chr --vcf {i} --make-bed "
                            f"--out {o} --pca 20 approx --update-sex {s} --max-alleles 2", shell=True)


if __name__ == '__main__':
    main(sys.argv[1:])
