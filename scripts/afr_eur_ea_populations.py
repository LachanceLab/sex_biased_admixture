#!/usr/bin/env python
import os
import sys
import argparse
import subprocess


def run_simulations(recipe, ots, rec, length, length_mtDNA):
    """
    Execute SLiM simulations
    :param recipe: str, path to recipe to use
    :param ots: str, path to output tree sequence
    :param rec: float, recombination rate
    :param length: int, chromosome length
    :param length_mtDNA: int, length of mtDNA
    """
    logfile = ots.replace('.trees', '.log')
    command = ' '.join(['slim', '-d', f"\"outfile='{ots}'\"", '-d', f"\"logfile='{logfile}'\"",
                        '-d', f'rec={rec}', '-d', f'L={length}', '-d', f'L_mtDNA={length_mtDNA}', recipe])
    subprocess.check_call(command, shell=True)


def main(argv):
    parser = argparse.ArgumentParser(description='Simulate ancestral continental mutations')
    parser.add_argument('-r', '--recombination_rate', type=float, default=1e-8, help='Recombination rate. [1e-8]')
    parser.add_argument('--recipe', help='SLiM recipe to use [scripts/afr_eur_ea.slim]', type=str,
                        default='scripts/afr_eur_ea.slim')
    parser.add_argument('-l', '--chromosome_length', help='Chromosome Length [1e6]', default=1e6, type=float)
    parser.add_argument('--mtDNA_length', help='Length of mitochondrial DNA [20000]', type=float, default=20000)
    parser.add_argument('-o', '--output_tree_sequences', help='File to write tree sequences to.'
                                                              '[./trees/afr_eur_ea.trees]',
                        default='./trees/afr_eur_ea.trees')
    args = parser.parse_args()
    output_dir = '/'.join(args.output_tree_sequences.split('/')[:-1])
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    run_simulations(args.recipe, args.output_tree_sequences, args.recombination_rate, args.chromosome_length,
                    args.mtDNA_length)


if __name__ == '__main__':
    main(sys.argv[1:])