#!/usr/bin/env python
import sys
import argparse
import subprocess


def run_simulations(recipe, ts, out, rec, length, mtDNA_length, afr, eur, ea, sm_afr, sm_eur, sm_ea, N_amr,
                    afr_const=None, eur_const=None, ea_const=None, sm_afr_const=None, sm_eur_const=None,
                    sm_ea_const=None, amr_growth_rate=0.0, rejection_rate_female_eur_male_non_eur=None):
    """
    Execute SLiM simulations
    :param recipe: str, path to SLiM recipe
    :param ts: str, path to input tree sequence from where to start simulations
    :param out: str, path where to write output tree sequence to
    :param rec: float, recombination rate
    :param length: int, chromosome length
    :param mtDNA_length: int, mtDNA length
    :param afr: float, African admixture fraction
    :param eur: float, European admixture fraction
    :param ea: float, East Asian admixture fraction
    :param sm_afr: float, fraction of African males
    :param sm_eur: float, fraction of European males
    :param sm_ea: float, fraction of East Asian males
    :param N_amr: float, population size of admixed population
    :param afr_const: float, African migration after initial admixture; default=None
    :param eur_const: float, European migration after initial admixture; default=None
    :param ea_const: float, East Asian migration after initial admixture; default=None
    :param sm_afr_const: float, fraction of constantly migrating African males; default=None
    :param sm_eur_const: float, fraction of constantly migrating European males; default=None
    :param sm_ea_const: float, fraction of constantly migrating East Asian males; default=None
    :param amr_growth_rate: float, exponential growth rate of admixed population
    :param rejection_rate_female_eur_male_non_eur: float, Fraction of matings of European females and non-European males
                                                          to reject in non_random mating scheme
    """
    logfile = f"{out.split('.trees')[0]}.log"
    command = ['slim', '-d', f"\"input_tree_sequence='{ts}'\"", '-d', f"\"outbase='{out}'\"",
               '-d', f"\"logfile='{logfile}'\"", '-d', f'rec={rec}', '-d', f'L={length}',
               '-d', f'L_mtDNA={mtDNA_length}', '-d', f'afr={afr}', '-d', f'eur={eur}', '-d', f'ea={ea}',
               '-d', f'sm_afr={sm_afr}',  '-d', f'sm_eur={sm_eur}',  '-d', f'sm_ea={sm_ea}',
               '-d', f'adx_n={N_amr}', '-d', f'amr_growth_rate={amr_growth_rate}']
    if not afr_const is None:
        command.extend(["-d", f'afr_const={afr_const}', "-d", f'eur_const={eur_const}', "-d", f'ea_const={ea_const}',
                        "-d", f'sm_afr_const={sm_afr_const}', "-d", f'sm_eur_const={sm_eur_const}',
                        "-d", f'sm_ea_const={sm_ea_const}'])
    if not rejection_rate_female_eur_male_non_eur is None and \
            (recipe == 'scripts/single_pulse_non_random_mating.slim' or
             recipe == 'scripts/constant_admixture_non_random_mating.slim'):
        command.extend(["-d", f'rejection_rate_female_eur_male_non_eur={rejection_rate_female_eur_male_non_eur}'])
    command.append(recipe)
    command = ' '.join(command)
    subprocess.check_call(command, shell=True)


def main(argv):
    parser = argparse.ArgumentParser(description="Simulate three-way (American) admixture of ancestral "
                                                 "continental populations")
    parser.add_argument('-r', '--recombination_rate', type=float, default=1e-8, help='Recombination rate. [1e-8]')
    parser.add_argument('--recipe', help='SLiM recipe to use [scripts/single_pulse_admixture.slim]', type=str,
                        default='scripts/single_pulse_admixture.slim')
    parser.add_argument('-l', '--chromosome_length', help='Chromosome Length [1e6]', default=1e6, type=float)
    parser.add_argument('--mtDNA_length', help='Length of mitochondrial DNA [20000]', type=float, default=20000)
    parser.add_argument('--tree_sequences', help='Input tree sequences [trees/afr_eur_ea.trees]',
                        default='trees/afr_eur_ea.trees')
    parser.add_argument('--afr', help='African admixture fraction [1/6]', default=1/6, type=float)
    parser.add_argument('--eur', help='European admixture fraction [1/3]', default=1/3, type=float)
    parser.add_argument('--ea', help='East Asian admixture fraction [1/2]', default=1/3, type=float)
    parser.add_argument('--sm_afr', help='Fraction of males of admixing Africans [1/2]', default=1/2, type=float)
    parser.add_argument('--sm_eur', help='Fraction of males of admixing Europeans[1/2]', default=1/2, type=float)
    parser.add_argument('--sm_ea', help='Fraction of males of admixing East Asians [1/2]', default=1/2, type=float)
    parser.add_argument('--N_amr', help='Initial size of admixed American population [1000]', default=1000, type=int)
    parser.add_argument('--afr_const', help='Constant African admixture fraction', default=None, type=float)
    parser.add_argument('--eur_const', help='Constant European admixture fraction', default=None, type=float)
    parser.add_argument('--ea_const', help='Constant East Asian admixture fraction', default=None, type=float)
    parser.add_argument('--sm_afr_const', help='Fraction of males of constantly admixing Africans', default=None,
                        type=float)
    parser.add_argument('--sm_eur_const', help='Fraction of males of constantly admixing Europeans', default=None,
                        type=float)
    parser.add_argument('--sm_ea_const', help='Fraction of males of constantly admixing East Asians', default=None,
                        type=float)
    parser.add_argument('--amr_growth_rate', help='Exponential growth rate of admixed populations [0.0]', default=0.0,
                        type=float)
    parser.add_argument('--rejection_rate_female_eur_male_non_eur',
                        help='Fraction of matings of European females and non-European males to reject in '
                             'non-random mating scheme. Only has an effect if '
                             '--recipe scripts/single_pulse_non_random_mating.slim', default=None, type=float)
    args = parser.parse_args()

    output_treeseq_base = f"{'/'.join(args.tree_sequences.split('/')[:-1])}/{args.recipe.split('/')[-1].split('.slim')[0]}_g"
    if not args.afr_const is None:
        assert args.eur_const is not None and args.ea_const is not None and args.sm_afr_const is not None and \
               args.sm_eur_const is not None and args.sm_ea_const is not None and "constant" in args.recipe, \
            "If specifing one constant admixture parameter you have to specify all and " \
            "select SLiM recipe for constant admixture"
    run_simulations(args.recipe, args.tree_sequences, output_treeseq_base, args.recombination_rate,
                    args.chromosome_length, args.mtDNA_length, args.afr, args.eur, args.ea, args.sm_afr,
                    args.sm_eur, args.sm_ea, args.N_amr, args.afr_const, args.eur_const, args.ea_const,
                    args.sm_afr_const, args.sm_eur_const, args.sm_ea_const, args.amr_growth_rate,
                    args.rejection_rate_female_eur_male_non_eur)


if __name__ == '__main__':
    main(sys.argv[1:])