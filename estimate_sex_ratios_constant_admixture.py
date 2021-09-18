#!/usr/bin/python
"""
#-----------------------------------------------------------
The code of the model has been adopted from Amy Goldberg's
Matlab code for the model (incl. comments)
https://www.genetics.org/content/201/1/263 (Eq. 19, 20 & 25)
-----------------------------------------------------------
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import argparse


def read_data(in_file):
    """
    Read data
    :param in_file: str, path to data file
    :return: pd.DataFrame, list, list, df contains X chromsomal and
                                       autosomal ancestry estimates for all populations and ancestries
    """
    df = pd.read_csv(in_file, sep=',', index_col=0, header=0, comment='#')
    populations = df.index
    ancestries = list(set([x.split('_')[0] for x in df.columns[:-1]]))
    return df, populations, ancestries


def constant_admixture(A, B, p, n=0.05):
    """
    Estimate female and male contribution under constant admixture.
    Assumes admixture occurred long enough ago (> 10 generations) so that limits have been reached
    :param A: float, estimated autosomal ancestry proportion
    :param B: float, estimated X chromosomal ancestry proportion
    :param p: float, fractions of females in sample
    :param n: float, increments for grid search
    :return: array, estimated parameters (s1f, s1m, s2f, s2m, hf, hm, euclidean distance)
    """
    # create grid for grid search
    s1f_vals, s1m_vals, s2f_vals, s2m_vals = np.meshgrid(np.arange(0.0, 1 + n, n),
                                                         np.arange(0.0, 1 + n, n),
                                                         np.arange(0.0, 1 + n, n),
                                                         np.arange(0.0, 1 + n, n))
    # avg contribution of source population 1
    s1 = (s1f_vals + s1m_vals) / 2
    # avg contribution of source population 2
    s2 = (s2f_vals + s2m_vals) / 2
    # female contribution of admixed pop
    hf_vals = np.ones_like(s1f_vals) - s1f_vals - s2f_vals
    # male contribution of admixed pop
    hm_vals = np.ones_like(s1m_vals) - s1m_vals - s2m_vals
    # avg contribution of admixed pop
    h = (hf_vals + hm_vals) / 2
    # conditions by Goldberg et al. --> admixed pop contributes to itself
    mask = np.where((hf_vals > 0.001) & (hm_vals > 0.001) & (h < 1))
    # select values that fulfill conditions
    s1 = s1[mask]
    s2 = s2[mask]
    h = h[mask]
    s1f_vals = s1f_vals[mask]
    s1m_vals = s1m_vals[mask]
    s2f_vals = s2f_vals[mask]
    s2m_vals = s2m_vals[mask]
    hf_vals = hf_vals[mask]
    hm_vals = hm_vals[mask]
    # expected autosomal ancestry proportions
    ExpH1 = s1 / (s1 + s2)
    # expected X chromosomal ancestry proportions in females (Eq. 19 in Goldberg et al.)
    ExpHf1X = (s1f_vals + s1m_vals + s1f_vals * hm_vals) /\
             (s1f_vals + s1m_vals + s2f_vals + s2m_vals + hm_vals * (s1f_vals + s2f_vals))
    # expected X chromosomal ancestry proportions in males (Eq 20 in Goldberg et al.)
    ExpHm1X = (2 * s1f_vals + hf_vals * s1m_vals) /\
             (s1f_vals + s1m_vals + s2f_vals + s2m_vals + hm_vals * (s1f_vals + s2f_vals))
    # calculates euclidean distance between observed data and calculations (Eq. 25 in Goldberg et al.)
    euc_distance = np.sqrt((ExpH1 - A) ** 2 + (p * ExpHf1X + (1 - p) * ExpHm1X - B) ** 2)
    # store in array
    output = np.stack([s1f_vals, s1m_vals, s2f_vals, s2m_vals, hf_vals, hm_vals, euc_distance]).T
    # sort by euclidean distance
    sorted = output[np.argsort(output[:, 6]), :]
    # percentile cutoff for parameters to keep, e.g. 1 is 1 percent
    m = 1
    i = round(m / 100 * output.shape[0])
    # accept top percentiles of parameter combinations
    accept = sorted[:i, :]
    return accept


def main(argv):
    parser = argparse.ArgumentParser(description="Re-implementations of Goldberg's recursion equation "
                                                 "to calculate the sex ratios for constant admixture "
                                                 "https://www.genetics.org/content/201/1/263 "
                                                 "Provided data file contains ancestry data from "
                                                 "Micheletti et al. (2020) "
                                                 "https://www.sciencedirect.com/science/article/pii/S0002929720302007?via%3Dihub#mmc6 "
                                                 "Sex ratios are inferred based on the provided data first, and then "
                                                 "for the corrected data. The corrected data distributes the missing"
                                                 "ancestries proportionally to all provided ancestries.")
    parser.add_argument('-i', '--input', help='Input file. First column contains index (target populations) and '
                                              'first row header (ancestries). Header should formatted so that acronym '
                                              'of ancestry comes first followed by "_" and either "A" or "X" depending '
                                              'on whether it refers to autosomal or X-related ancestry. For instance, '
                                              'for african ancestry the header should contain the columns AFR_X '
                                              'and AFR_A. For each ancestry, autosomal and X-related ancestry must be '
                                              'given. The last column should be named "pf" and represent proportion of '
                                              'females contributing to the gene pool. Comments can be included in '
                                              'separate lines, starting with a "#", and are ignored while parsing. '
                                              'The file should be comma separated.')
    parser.add_argument('-o', '--output', help='Output filename. Will be an excel file.')
    parser.add_argument('-n', '--increment_grid_search', help='Increment for gird search to find parameters with lowest'
                                                              'Euclidean distance, default=0.05', default=0.05,
                        type=float)
    parser.add_argument('-p', '--proportionally', help='Set flag to distribute missing ancestries proportionally.'
                                                       'Default is equally', action='store_true', default=False,
                        required=False)
    args = parser.parse_args()
    df, populations, ancestries = read_data(args.input)
    # raw values
    fig, ax = plt.subplots(len(populations), 1, sharex=True, figsize=(10, 10))
    plt.subplots_adjust(hspace=0.65)
    # accoutned for missing ancestries
    fig1, ax1 = plt.subplots(len(populations), 1, sharex=True, figsize=(10, 10))
    plt.subplots_adjust(hspace=0.65)

    for n, population in enumerate(populations):
        for i, ancestry in enumerate(ancestries):
            parameter_estimates = constant_admixture(df.loc[population, '{}_A'.format(ancestry)],
                                                     df.loc[population, '{}_X'.format(ancestry)],
                                                     df.loc[population, 'pf'], args.increment_grid_search)
            boxplot = ax[n].boxplot(np.log(parameter_estimates[:, 0] / parameter_estimates[:, 1]), positions=[i],
                                    vert=False, widths=0.5, notch=True, patch_artist=True,
                                    whis=[2.5, 97.5], showfliers=False)
            for patch, color in zip(boxplot['boxes'], ['red']):
                patch.set_facecolor(color)
        ax[n].set_yticks(range(0, len(ancestries)))
        ax[n].set_yticklabels(ancestries)
        ax[n].set_title(population)
        # correct data --> distribute missing ancestry proportionally to all provided ancestries
        if args.proportionally:
            x_contributions = np.array([df.loc[population, '{}_X'.format(a)] for a in ancestries])
            a_contributions = np.array([df.loc[population, '{}_A'.format(a)] for a in ancestries])
            x_contributions_adj = x_contributions / x_contributions.sum()
            a_contributions_adj = a_contributions / a_contributions.sum()
        # correct data --> distribute missing ancestry equally to all provided ancestries
        else:
            x_contributions = np.array([df.loc[population, '{}_X'.format(a)] for a in ancestries])
            a_contributions = np.array([df.loc[population, '{}_A'.format(a)] for a in ancestries])
            x_contributions_adj = x_contributions + (1 - x_contributions.sum()) / x_contributions.shape[0]
            a_contributions_adj = a_contributions + (1 - a_contributions.sum()) / a_contributions.shape[0]

        for i, (x, a, ancestry) in enumerate(zip(x_contributions_adj, a_contributions_adj, ancestries)):
            parameter_estimates = constant_admixture(a, x, df.loc[population, 'pf'], args.increment_grid_search)
            boxplot = ax1[n].boxplot(np.log(parameter_estimates[:, 0] / parameter_estimates[:, 1]), positions=[i],
                                     vert=False, widths=0.5, notch=True, patch_artist=True,
                                     whis=[2.5, 97.5], showfliers=False)
            for patch, color in zip(boxplot['boxes'], ['red']):
                patch.set_facecolor(color)
        ax1[n].set_yticks(range(0, len(ancestries)))
        ax1[n].set_yticklabels(ancestries)
        ax1[n].set_title(population)
    ax[n].set_xlabel(r"$ln(s_f/s_m)$")
    ax[n].set_xlim([-4, 4])
    fig.savefig('sex_specific_contribution_constant_admixture.jpeg', bbox_inches='tight')
    ax1[n].set_xlabel(r"$ln(s_f/s_m)$")
    ax1[n].set_xlim([-4, 4])
    if args.proportionally:
        fig1.savefig('sex_specific_contribution_constant_admixture_proportionally.jpeg', bbox_inches='tight')
    else:
        fig1.savefig('sex_specific_contribution_constant_admixture_equally.jpeg', bbox_inches='tight')


if __name__ == '__main__':
    main(sys.argv[1:])
