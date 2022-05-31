#!/usr/bin/env python
"""
---------------------------------------------------------
Re-implementations of Goldberg's recursion equation
to calculate the sex ratios for a single pulse admixture
https://www.genetics.org/content/201/1/263
---------------------------------------------------------
"""
import argparse
import sys
import numpy as np
import pandas as pd
from sklearn.utils import resample


def read_results_admixture(admixture_results_autosomes, admixture_results_chrx, fam_file):
    """
    Read Admixture results for autosomes and X chromosomes as well fam file to obtain individual IDs and sexes
    :param admixture_results_autosomes: str, Path to *.Q file for autosomes
    :param admixture_results_chrx: str, Path to *.Q file for X chromosomes
    :param fam_file: str, Path to *.fam file
    :return: pd.DataFrame, pd.DataFrame, admixture results for autosomes and X chromosomes
    """
    # read fam file
    indv_names_and_sex = pd.read_csv(fam_file, sep='\t', header=None, usecols=[1, 4])
    # read Q file autosomes
    K = int(admixture_results_autosomes.split('.')[-2])
    admixture_autosomes = pd.read_csv(admixture_results_autosomes, sep=' ', header=None,
                                      names=[f'K{i}' for i in range(K)])
    # send IDs as index
    admixture_autosomes.index = indv_names_and_sex.iloc[:, 0].values
    # add sexes
    admixture_autosomes.loc[:, 'sex'] = indv_names_and_sex.iloc[:, 1].values - 1
    # read Q file X chromosomes
    admixture_x = pd.read_csv(admixture_results_chrx, sep=' ', header=None, names=[f'K{i}' for i in range(K)])
    # set IDs as index
    admixture_x.index = indv_names_and_sex.iloc[:, 0].values
    # add sexes
    admixture_x.loc[:, 'sex'] = indv_names_and_sex.iloc[:, 1].values - 1
    return admixture_autosomes, admixture_x


def read_data_summary_statistics(in_file):
    """
    Read data
    :param in_file: str, path to data file
    :return: pd.DataFrame, list, list, df contains X chromsomal and
                                       autosomal ancestry estimates for all populations and ancestries
    """
    df = pd.read_csv(in_file, sep=',', index_col=0, header=0, comment='#')
    populations = df.index
    ancestries = list(sorted(set([x.split('_')[0] for x in df.columns[:-1]])))
    return df, populations, ancestries


def calculate_sex_contributions(x_contribution, autosomal_contribution, generations, pf):
    """
    Calculate sex ratios for specified generations
    :param x_contribution: float, inferred x chromosomal ancestry proportions
    :param autosomal_contribution: float, inferred autosomal ancestry proportions
    :param generations: array-like, generations (must be greater than or equal to 2)
    :param pf: float, proportions of females in sample
    :return: (array-like, array-like, float, float) female contribution per generation,
                                                    male contribution per generation
                                                    female contribution in equilibrium
                                                    male contribution in equilibrium
    """
    # calculate sf and sm in the limit of infinity
    admixture = np.array([x_contribution, autosomal_contribution])
    transitions = np.array([[2 / 3, 1 / 3],
                            [1 / 2, 1 / 2]])
    sf_inf, sm_inf = np.linalg.solve(transitions, admixture)
    # calculate sf and sm for more recent admixture
    s_f = []
    s_m = []
    for i, g in enumerate(generations):
        # weigh by female contribution (pf) and male contribution (1 - pf), pf is the proportion of females in sample
        # X-related ancestry according to equation 12 and 13 in Goldberg
        # autosomal ancestry is just 1/2 sf + 1/2 sm
        # [[contribution of females to X ancestry in females + contribution of females to X ancestry in males,
        #   contribution of males to X ancestry in females + contribution of males to X ancestry in males],
        #  [contribution of females to autosomal ancestry, contribution of males to autosomal ancestry]]
        transitions = np.array([[((2 + (-0.5) ** g) / 3 * pf + (2 + (-0.5) ** (g - 1)) / 3 * (1 - pf)),
                                 ((1 - (-0.5) ** g) / 3 * pf + (1 - (-0.5) ** (g - 1)) / 3 * (1 - pf))],
                                [1 / 2, 1 / 2]])
        # solve system of linear equations for sf and sm
        c_s_f, c_s_m = np.linalg.solve(transitions, admixture)
        s_f.append(c_s_f)
        s_m.append(c_s_m)

    return np.array(s_f), np.array(s_m), sf_inf, sm_inf


def compute_sex_ratios_from_summary_statistics(input_file, output_file):
    """
    Compute sex ratios from summary statistics (mean ancestry proportions)
    :param input_file: str, Path to file with mean ancestry proportions
    :param output_file: str, Path where to write results to
    """
    df, populations, ancestries = read_data_summary_statistics(input_file)

    generations = np.arange(2, 16, 1).tolist()
    generations.append('inf')
    complete_df = pd.DataFrame()
    i = 0
    for population in populations:
        c_df = pd.DataFrame()
        for ancestry in ancestries:
            # calculate sex ratios using Goldberg approach using input data
            sf, sm, sf_inf, sm_inf = calculate_sex_contributions(df.loc[population, '{}_X'.format(ancestry)],
                                                                 df.loc[population, '{}_A'.format(ancestry)],
                                                                 generations[:-1],
                                                                 df.loc[population, 'pf'])
            sf = np.append(sf, sf_inf)
            sm = np.append(sm, sm_inf)
            c_df['population'] = [population] * len(sf)
            c_df['G'] = generations
            c_df['{}_X'.format(ancestry)] = [df.loc[population, '{}_X'.format(ancestry)]] * len(sf)
            c_df['{}_A'.format(ancestry)] = [df.loc[population, '{}_A'.format(ancestry)]] * len(sf)
            # do rounding of sex-specific contributions that are added to df --> 3 decimals
            c_df['{}_sf'.format(ancestry)] = sf.round(3)
            c_df['{}_sm'.format(ancestry)] = sm.round(3)
            # calculate sex ratios with exact numbers
            females_to_male = sf / sm
            males_to_female = sm / sf
            # round sex ratios to three significant digits
            c_df['{}_sf/sm'.format(ancestry)] = [x.round(2) if np.abs(x) < 10 else x.round(1) if np.abs(x) > 10
                                                                                                 and np.abs(
                x) < 100 else int(x.round(0)) for x in females_to_male]
            c_df['{}_sm/sf'.format(ancestry)] = [x.round(2) if np.abs(x) < 10 else x.round(1) if np.abs(x) > 10
                                                                                                 and np.abs(
                x) < 100 else int(x.round(0)) for x in males_to_female]
        complete_df = complete_df.append(c_df)
    # save to excel
    complete_df.reset_index(drop=True, inplace=True)
    complete_df.to_excel(output_file, index=False)


def bootstrap_raw_proportions(admixture_results_autosomes, admixture_results_chrx, fam_file, output_file,
                              n_resamples=None):
    """
    Bootstrap raw ancestry estimates from ADMIXTURE output files to obtain confidence intervals of mean
    ancestry estimates and sex ratios
    :param admixture_results_autosomes: str, Path to .Q file for autosomes
    :param admixture_results_chrx: str, Path to .Q file for X chromosomes
    :param fam_file: str, path to corresponding .fam file
    :param output_file: str, Path where to write results to
    :param n_resamples: int, number of resampling steps
    """
    # generations for Goldberg model
    generations = np.arange(2, 16, 1).tolist()
    # equilibrium model
    generations.append('inf')
    # read raw data
    admixture_autosomes, admixture_x = read_results_admixture(admixture_results_autosomes, admixture_results_chrx,
                                                              fam_file)
    # join
    ancestry_proportions = admixture_autosomes.join(admixture_x.drop(["sex"], axis=1), lsuffix='_A',
                                                    rsuffix='_X')
    # get ancestries
    ancestries = [col.split('_')[0] for col in ancestry_proportions.columns if col.endswith('A')]
    # get sample names
    sample_names = ancestry_proportions.index.values
    # bootstrap
    if n_resamples is None:
        n_resamples = max([min([ancestry_proportions.shape[0], 10000]), 500])
    iterations_df = pd.DataFrame(index=np.repeat(np.arange(0, n_resamples), len(generations)))

    for i in range(n_resamples):
        # resample
        resampled = resample(sample_names)
        resampled_autosomes = ancestry_proportions.loc[resampled, [f"{anc}_A" for anc in ancestries]]
        resampled_x = ancestry_proportions.loc[resampled, [f"{anc}_X" for anc in ancestries]]
        mean_sex = ancestry_proportions.loc[resampled, "sex"].mean()
        # iterate over ancestries
        for ancestry in ancestries:
            # calculate mean proportions
            mean_autosomal_prop = resampled_autosomes.loc[:, f"{ancestry}_A"].mean()
            mean_x_prop = resampled_x.loc[:, f"{ancestry}_X"].mean()
            # get sex specific contributions
            sf, sm, sf_inf, sm_inf = calculate_sex_contributions(mean_x_prop, mean_autosomal_prop, generations[:-1],
                                                                 mean_sex)
            sf = np.append(sf, sf_inf)
            sm = np.append(sm, sm_inf)
            iterations_df.loc[i, 'G'] = generations
            iterations_df.loc[i, '{}_X'.format(ancestry)] = [mean_x_prop] * len(sf)
            iterations_df.loc[i, '{}_A'.format(ancestry)] = [mean_autosomal_prop] * len(sf)
            # do rounding of sex-specific contributions that are added to df --> 3 decimals
            iterations_df.loc[i, '{}_sf'.format(ancestry)] = sf.round(3)
            iterations_df.loc[i, '{}_sm'.format(ancestry)] = sm.round(3)
            # calculate sex ratios with exact numbers
            females_to_male = sf / sm
            males_to_female = sm / sf
            # round sex ratios to three significant digits
            iterations_df.loc[i, '{}_sf/sm'.format(ancestry)] = [x.round(2) if np.abs(x) < 10
                                                                 else x.round(1) if np.abs(x) > 10 and np.abs(x) < 100
                                                                 else int(x.round(0)) for x in females_to_male]
            iterations_df.loc[i, '{}_sm/sf'.format(ancestry)] = [x.round(2) if np.abs(x) < 10
                                                                 else x.round(1) if np.abs(x) > 10 and np.abs(x) < 100
                                                                 else int(x.round(0)) for x in males_to_female]
    # summarize
    mean_vals = iterations_df.groupby("G").mean()
    # 95% CI
    ci_low = iterations_df.groupby('G').quantile(0.025)
    ci_up = iterations_df.groupby('G').quantile(0.975)
    ci_up.columns = [f"{col}_ci_up" for col in ci_up.columns]
    summary_df = ci_low.join(mean_vals, lsuffix="_ci_low", rsuffix="_mean").join(ci_up)
    # save
    summary_df = summary_df.loc[:, np.sort(summary_df.columns)]
    summary_df.to_excel(output_file, index=True)


def main(argv):
    parser = argparse.ArgumentParser(description="Re-implementations of Goldberg's recursion equation "
                                                 "to calculate the sex ratios for a single pulse admixture "
                                                 "https://www.genetics.org/content/201/1/263 "
                                                 "Provided data file contains ancestry estimates from"
                                                 "Micheletti et al. (2020)"
                                                 "https://doi.org/10.1016/j.ajhg.2020.06.012 "
                                                 "If mean ancestry proportions are provided (-i), sex ratios are "
                                                 "inferred based on the provided data first,"
                                                 "and then after distributing unassigned ancestries in all possible"
                                                 " ways such that the ancestries sum up to 100%. "
                                                 "If .Q files are provided the individuals are bootstrapped, allowing "
                                                 "to report confidence intervals for sex specific contributions")
    parser.add_argument('-i', '--input', help='Input file with mean ancestry proportions. '
                                              'First column contains index (target populations) and '
                                              'first row header (ancestries). Header should formatted so that acronym '
                                              'of ancestry comes first followed by "_" and either "A" or "X" depending '
                                              'on whether it refers to autosomal or X-related ancestry. For instance, '
                                              'for african ancestry the header should contain the columns AFR_X '
                                              'and AFR_A. For each ancestry, autosomal and X-related ancestry must be '
                                              'given. The last column should be named "pf" and represent proportion of '
                                              'females contributing to the gene pool. Comments can be included in '
                                              'separate lines, starting with a "#", and are ignored while parsing. '
                                              'The file should be comma separated.', required=False)
    parser.add_argument('-qa', '--admixture_results_autosomes', help='.Q file from Admixture for autosomes',
                        required=False)
    parser.add_argument('-qx', '--admixture_results_chrx', help='.Q file from Admixture for X chromosome',
                        required=False)
    parser.add_argument('-f', '--fam_file', help='*.fam file generated by plink to get individual IDs', required=False)
    parser.add_argument('-b', '--bootstraps', default=None,
                        help='Number of bootstrap resamples. Default is the number of samples '
                             'but at least 500 and at most 10000', type=int)
    parser.add_argument('-o', '--output', help='Output filename. Will be an excel file.', required=True)
    args = parser.parse_args()
    if args.input:
        compute_sex_ratios_from_summary_statistics(args.input, args.output)
    elif args.admixture_results_autosomes and args.admixture_results_chrx and args.fam_file:
        bootstrap_raw_proportions(args.admixture_results_autosomes, args.admixture_results_chrx, args.fam_file,
                                  args.output, args.bootstraps)
    else:
        raise AttributeError("Either specify summary statistic file with -i flag or Admixture *Q file for autosomomes "
                             "and X chromsomes and fam file using -qa, -qx and -f flag, respectively.")


if __name__ == '__main__':
    main(sys.argv[1:])
