#!/usr/bin/python
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
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Polygon


def read_data(in_file):
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


def sensitivity_analysis(autosomal_ancestry):
    """
    Plot sensitivity analysis to show divergence (Figure 1)
    :param autosomal_ancestry: float, autosomal ancestry proportion
    """
    # function for calculating equilibrium sex ratio
    sf_sm = lambda x, a: (3 * x - 2 * a) / (4 * a - 3 * x)
    x_ancestry = np.arange(0, 1.001, 0.001)
    # compute sex ratios for different x ancestry proportions
    sex_ratio = sf_sm(x_ancestry, autosomal_ancestry)
    sex_ratio[sex_ratio > 500] = np.inf
    sex_ratio[sex_ratio < -500] = -np.inf
    fig, ax = plt.subplots(1, 2, figsize=(7.5, 3.75))
    plt.subplots_adjust(wspace=0.35)
    ax[0].xaxis.set_minor_locator(MultipleLocator(0.1))
    ax[0].xaxis.set_major_locator(MultipleLocator(0.2))
    ax[0].yaxis.set_major_locator(MultipleLocator(2))
    ax[0].yaxis.set_minor_locator(MultipleLocator(1))
    # make plot squared
    ax[0].set_aspect(0.05)
    ax[0].plot(x_ancestry, sex_ratio, lw=1)
    ax[0].set_xlabel(r'X chromosomal ancestry proportion $\left(H^X_1\right)$', labelpad=10, fontdict=dict(fontsize=7))
    ax[0].set_ylabel(r"Ratio of females to one male $\left(s^f_1/s^m_1\right)$", fontdict=dict(fontsize=7))
    ax[0].set_ylim([-10, 10])
    ax[0].set_xlim([0, 1.0])
    ax[0].annotate(r'$H^A_1$', xy=(autosomal_ancestry, -0.005), xytext=(autosomal_ancestry, -0.07),
                   arrowprops=dict(width=0.2, headwidth=1.75, headlength=1.75, facecolor='black'),
                   xycoords='axes fraction', fontsize=6, horizontalalignment="center")
    ax[0].axhline(1, c='black', ls='--', label=r'$s^f_1/s^m_1=1$', lw=0.75)
    ax[0].axhline(0, c='grey', ls=':', label=r'$s^f_1/s^m_1=0$', lw=0.75)
    ax[0].tick_params(labelsize=6)
    ax[0].legend(fontsize=7)
    ax[0].text(-0.22, ax[0].get_ylim()[1] + ax[0].get_ylim()[1] * 0.08, 'A', fontsize=12, weight='bold')
    # compute patches for different regimes
    a_vals = np.arange(0.00000001, 2.00000001, 0.00000001)
    sex_ratios_low_autosomal_prop = sf_sm(0.0000001, a_vals)
    sex_ratios_high_autosomal_prop = sf_sm(0.9999999, a_vals)

    # negative ratios due to too little X
    patch_neg_too_little_x = [[0.0000001, a_vals[np.where((sex_ratios_low_autosomal_prop < 0) &
                                                          (a_vals > 0.0000001))[0][0]]],  # bottom left
                              [0.99999999, a_vals[np.where((sex_ratios_high_autosomal_prop < 0) &
                                                           (a_vals > 0.9999999))[0][0]]],  # top right
                              [0, 1]]  # top left
    # male-biased admixture
    patch_male = [[0.0000001, a_vals[np.where((sex_ratios_low_autosomal_prop < 1) &
                                              (sex_ratios_low_autosomal_prop > 0) &
                                              (a_vals > 0.0000001))[0][0]]],  # bottom left
                  patch_neg_too_little_x[1],  # top left
                  [0.99999999, a_vals[np.where((sex_ratios_high_autosomal_prop < 1) &
                                               (a_vals > 0.9999999) &
                                               (sex_ratios_high_autosomal_prop > 0))[0][0]]]  # top right
                  ]
    # female-biased admixture
    patch_female = [[0.0000001, a_vals[np.where((sex_ratios_low_autosomal_prop > 1) &
                                                (a_vals < 0.0000001))[0][0]]],  # bottom left
                    [0.99999999, a_vals[np.where((sex_ratios_high_autosomal_prop > 1) &
                                                 (a_vals < 0.9999999))[0][0]]],  # top right
                    patch_male[1]]  # top left

    # negative sex ratios due to too much X
    patch_neg_too_much_x = [[0.0000001, a_vals[np.where((sex_ratios_low_autosomal_prop < 0) &
                                                        (a_vals < 0.0000001))[0][0]]],  # bottom left
                            patch_female[1],  # top
                            [1, 0]]  # bottom right
    # plot patches
    ax[1].xaxis.set_major_locator(MultipleLocator(0.2))
    ax[1].xaxis.set_minor_locator(MultipleLocator(0.1))
    ax[1].yaxis.set_major_locator(MultipleLocator(0.2))
    ax[1].yaxis.set_minor_locator(MultipleLocator(0.1))
    ax[1].add_patch(Polygon(patch_neg_too_much_x, closed=True, fill=True, color='darkgray', label='model breakage'))
    ax[1].add_patch(Polygon(patch_female, closed=True, fill=True, color='cornflowerblue',
                            label='female-biased'))
    ax[1].add_patch(Polygon(patch_male, closed=True, fill=True, color='lightskyblue', label='male-biased'))
    ax[1].add_patch(Polygon(patch_neg_too_little_x, closed=True, fill=True, color='darkgray'))
    ax[1].plot([0, 1], [0, 1], ls='--', color='black', lw= 0.75)
    ax[1].set_xlim([0, 1.0])
    ax[1].set_ylim([0, 1.0])
    ax[1].set_xlabel(r'X chromosomal ancestry proportion $\left(H^X_1\right)$', labelpad=10, fontdict=dict(fontsize=7))
    ax[1].set_ylabel(r'Autosomal ancestry proportion $\left(H^A_1\right)$', fontdict=dict(fontsize=7))
    ax[1].tick_params(labelsize=6)
    ax[1].set_aspect('equal')

    ax[1].legend(bbox_to_anchor=(1.025, 0.5), loc='center left', fontsize=7)
    ax[1].text(-0.22, ax[1].get_ylim()[1] + ax[1].get_ylim()[1] * 0.04, 'B', fontsize=12, weight='bold')

    fig.savefig('figure1.pdf', dpi=500, bbox_inches='tight')
    plt.close()


def distribute_ancestries_proportionally(df, population, ancestries):
    """
    Distribute the unassigned ancestries proportionally into all bins
    :param df: pd.DataFrame, X chromosomal and autosomal ancestry estimates
    :param population: str, current population
    :param ancestries: list, ancestries to assign missing ancestries to
    :return: tuple, (adjusted X chromosomal ancestry estimates, adjusted autosomal ancestry estimates)
                    in order of ancestries
    """
    # correct data --> distribute missing ancestry proportionally to all provided ancestries
    x_contributions = np.array([df.loc[population, '{}_X'.format(a)] for a in ancestries])
    a_contributions = np.array([df.loc[population, '{}_A'.format(a)] for a in ancestries])
    x_contributions_adj = x_contributions / x_contributions.sum()
    a_contributions_adj = a_contributions / a_contributions.sum()
    return np.round(x_contributions_adj, 4), np.round(a_contributions_adj, 4)


def distribute_ancestries_equally(df, population, ancestries):
    """
    Distribute the unassigned ancestries equally into all bins
    :param df: pd.DataFrame, X chromosomal and autosomal ancestry estimates
    :param population: str, current population
    :param ancestries: list, ancestries to assign missing ancestries to
    :return: tuple, (adjusted X chromosomal ancestry estimates, adjusted autosomal ancestry estimates)
                    in order of ancestries
    """
    # correct data --> distribute missing ancestry equally to all provided ancestries
    x_contributions = np.array([df.loc[population, '{}_X'.format(a)] for a in ancestries])
    a_contributions = np.array([df.loc[population, '{}_A'.format(a)] for a in ancestries])
    x_contributions_adj = x_contributions + (1 - x_contributions.sum()) / x_contributions.shape[0]
    a_contributions_adj = a_contributions + (1 - a_contributions.sum()) / a_contributions.shape[0]
    return np.round(x_contributions_adj, 4), np.round(a_contributions_adj, 4)


def main(argv):
    parser = argparse.ArgumentParser(description="Re-implementations of Goldberg's recursion equation "
                                                 "to calculate the sex ratios for a single pulse admixture "
                                                 "https://www.genetics.org/content/201/1/263 "
                                                 "Provided data file contains ancestry estimates from"
                                                 "Micheletti et al. (2020)"
                                                 "https://www.sciencedirect.com/science/article/pii/S0002929720302007?via%3Dihub#mmc6 "
                                                 "Sex ratios are inferred based on the provided data first,"
                                                 "and then for the corrected data. The corrected data distributes "
                                                 "the missing ancestries equally or proportionally to all "
                                                 "provided ancestries.")
    parser.add_argument('-i', '--input', help='Input file. First column contains index (target populations) and '
                                              'first row header (ancestries). Header should formatted so that acronym '
                                              'of ancestry comes first followed by "_" and either "A" or "X" depending '
                                              'on whether it refers to autosomal or X-related ancestry. For instance, '
                                              'for african ancestry the header should contain the columns AFR_X '
                                              'and AFR_A. For each ancestry, autosomal and X-related ancestry must be '
                                              'given. The last column should be named "pf" and represent proportion of '
                                              'females contributing to the gene pool. Comments can be included in '
                                              'separate lines, starting with a "#", and are ignored while parsing. '
                                              'The file should be comma separated.', required=True)
    parser.add_argument('-o', '--output', help='Output filename. Will be an excel file.', required=True)
    parser.add_argument('-p', '--proportionally', help='Set flag to distribute missing ancestries proportionally.'
                                                       'Default is equally', action='store_true', default=False,
                        required=False)
    args = parser.parse_args()
    df, populations, ancestries = read_data(args.input)
    sensitivity_analysis(0.123)
    generations = np.arange(2, 16, 1).tolist()
    generations.append('inf')
    complete_df = pd.DataFrame()
    complete_df_adj = pd.DataFrame()

    for population in populations:
        c_df = pd.DataFrame()
        c_df_adj = pd.DataFrame()
        for ancestry in ancestries:
            # calculate sex ratios using Goldberg approach using input data
            sf, sm, sf_inf, sm_inf = calculate_sex_contributions(df.loc[population, '{}_X'.format(ancestry)],
                                                                 df.loc[population, '{}_A'.format(ancestry)], generations[:-1],
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
                                                 and np.abs(x) < 100 else int(x.round(0)) for x in females_to_male]
            c_df['{}_sm/sf'.format(ancestry)] = [x.round(2) if np.abs(x) < 10 else x.round(1) if np.abs(x) > 10
                                                         and np.abs(x) < 100 else int(x.round(0)) for x in males_to_female]
        if args.proportionally:
            x_contributions_adj, a_contributions_adj = distribute_ancestries_proportionally(df, population, ancestries)
        else:
            x_contributions_adj, a_contributions_adj = distribute_ancestries_equally(df, population, ancestries)
        for x, a, ancestry in zip(x_contributions_adj, a_contributions_adj, ancestries):
            # calculate sex ratios using Goldberg approach using corrected data
            # missing ancestry is distributed proportional to all provided ancestries.
            sf, sm, sf_inf, sm_inf = calculate_sex_contributions(x, a, generations[:-1], df.loc[population, 'pf'])
            sf = np.append(sf, sf_inf)
            sm = np.append(sm, sm_inf)
            c_df_adj['population'] = [population] * len(sf)
            c_df_adj['G'] = generations
            c_df_adj['{}_X_adj'.format(ancestry)] = [x] * len(sf)
            c_df_adj['{}_A_adj'.format(ancestry)] = [a] * len(sf)
            # do rounding of sex-specific contributions  that are added to df --> 3 decimals
            c_df_adj['{}_sf_adj'.format(ancestry)] = sf.round(3)
            c_df_adj['{}_sm_adj'.format(ancestry)] = sm.round(3)
            # calculate sex ratios with exact numbers
            females_to_male = sf / sm
            males_to_female = sm / sf
            # round sex ratios to three significant digits
            c_df_adj['{}_sf/sm_adj'.format(ancestry)] = [x.round(2) if np.abs(x) < 10 else x.round(1) if np.abs(x) > 10
                                                         and np.abs(x) < 100 else int(x.round(0)) for x in females_to_male]
            c_df_adj['{}_sm/sf_adj'.format(ancestry)] = [x.round(2) if np.abs(x) < 10 else x.round(1) if np.abs(x) > 10
                                                         and np.abs(x) < 100 else int(x.round(0)) for x in males_to_female]
        complete_df = complete_df.append(c_df)
        complete_df_adj = complete_df_adj.append(c_df_adj)

    # save to excel
    complete_df.reset_index(drop=True, inplace=True)
    complete_df.to_excel(args.output, index=False)
    complete_df_adj.reset_index(drop=True, inplace=True)
    if args.proportionally:
        complete_df_adj.to_excel(args.output.replace('.xlsx', '_adjusted_proportionally.xlsx'), index=False)
    else:
        complete_df_adj.to_excel(args.output.replace('.xlsx', '_adjusted_equally.xlsx'), index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
