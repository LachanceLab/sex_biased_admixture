#!/usr/bin/env python
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


def expected_x_proportions(s1f_0, s1m_0, s1f, s1m, s2f_0, s2m_0, s2f, s2m, target_g=15):
    """
    Compute expectation of X chromosomal ancestry proportion in females and males using Equations A2 and A3 in
    Goldberg and Rosenberg (2014)
    :param s1f_0: float, initial female contribution from source population 1
    :param s1m_0: float, initial male contribution from source population 1
    :param s1f: float, constant female contribution from source population 1
    :param s1m: float, constant male contribution from source population 1
    :param s2f_0: float, initial female contribution from source population 2
    :param s2m_0: float, initial male contribution from source population 2
    :param s2f: float, constant female contribution from source population 2
    :param s2m: float, constant male contribution from source population 2
    :param target_g: int, generations after admixture
    :return: (float, float), expected X chromosomal ancestry proportion in females and males
    """
    # compute X chromosome proportions in generation 1
    hxm = s1f_0 / (s1f_0 + s1m_0 + s2f_0 + s2m_0)
    hxf = (s1f_0 + s1m_0) / (s1f_0 + s1m_0 + s2f_0 + s2m_0)
    hxm_vals = [hxm]
    hxf_vals = [hxf]
    # compute sex specific contributions from admixed population
    hf = 1 - (s1f + s2f)
    hm = 1 - (s1m + s2m)
    # recursively compute ancestry proportions until target generation is reached
    for g in range(1, target_g):
        new_hxm = s1f + hf * hxf_vals[g - 1]
        if g == 1:
            new_hxf = (s1f + s1m) / 2 + 0.5 * (hxf_vals[g - 1] * hf + s1f_0 * hm)
        else:
            new_hxf = (s1f + s1m) / 2 + 0.5 * hf * hxf_vals[g - 1] + (0.5 * hm) * (s1f + hf * hxf_vals[g - 2])
        hxf_vals.append(new_hxf)
        hxm_vals.append(new_hxm)
    # return last generation
    return hxf_vals[-1], hxm_vals[-1]


def expected_autosome_proportions(s1_0, s1, h, g):
    """
    Compute the expected autosomal ancestry proportion. Equation 37 in Goldberg, Verdu and Rosenberg (2014)
    :param s1_0: float, initial contribution form source population 1
    :param s1: float, constant contribution from source population 1
    :param h: float, contribution from admixed population
    :param g: int, generation after admixture
    :return: float, autosomal admixture proportion
    """
    return s1_0 * h ** (g - 1) + s1 * ((1 - h ** (g - 1)) / (1 - h))


def constant_admixture(A, B, p, n=0.02, threshold=0.01, g=15, s1f_0=1/2, s1m_0=1/2):
    """
    Estimate female and male contribution under constant admixture.
    Assumes admixture occurred long enough ago (> 10 generations) so that limits have been reached
    :param A: float, estimated autosomal ancestry proportion
    :param B: float, estimated X chromosomal ancestry proportion
    :param p: float, fractions of females in sample
    :param n: float, increments for grid search
    :param threshold: float, max Euclidean distance to accept parameter combination
    :param g: int, generations after admixture
    :param s1f_0: float, initial female contribution from source population 1
    :param s1m_0: float, initial male contribution from source population 1
    :return: array, estimated parameters (s1f, s1m, s2f, s2m, hf, hm, euclidean distance)
    """
    # create grid for grid search
    s1f_vals, s1m_vals, s2f_vals, s2m_vals = np.meshgrid(np.arange(0.0, 1 + n, n),
                                                         np.arange(0.0, 1 + n, n),
                                                         np.arange(0.0, 1 + n, n),
                                                         np.arange(0.0, 1 + n, n))
    # avg contribution of source population 1
    s1 = (s1f_vals + s1m_vals) / 2
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
    h = h[mask]
    s1f_vals = s1f_vals[mask]
    s1m_vals = s1m_vals[mask]
    s2f_vals = s2f_vals[mask]
    s2m_vals = s2m_vals[mask]
    hf_vals = hf_vals[mask]
    hm_vals = hm_vals[mask]
    # expected autosomal ancestry proportions
    ExpH1 = expected_autosome_proportions((s1f_0 + s1m_0) / 2, s1, h, g)
    # expected X chromosomal ancestry proportions in females (Eq. 19 in Goldberg et al.)
    ExpHf1X, ExpHm1X = expected_x_proportions(s1f_0, s1m_0, s1f_vals, s1m_vals, 1 - s1f_0, 1 - s1m_0,
                                              s2f_vals, s2m_vals, g)
    # calculates euclidean distance between observed data and calculations (Eq. 25 in Goldberg et al.)
    euc_distance = np.sqrt((ExpH1 - A) ** 2 + (p * ExpHf1X + (1 - p) * ExpHm1X - B) ** 2)
    # store in array
    output = np.stack([s1f_vals, s1m_vals, s2f_vals, s2m_vals, hf_vals, hm_vals, euc_distance]).T
    # cut-off by euclidean distance
    accept = output[output[:, 6] <= threshold]
    return accept


def main(argv):
    parser = argparse.ArgumentParser(description="Re-implementations of Goldberg's recursion equation "
                                                 "to calculate the sex ratios for constant admixture "
                                                 "https://www.genetics.org/content/201/1/263 "
                                                 "Provided data file contains ancestry data from "
                                                 "Micheletti et al. (2020) "
                                                 "https://doi.org/10.1016/j.ajhg.2020.06.012 "
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
                                              'The file should be comma separated.', required=True)
    parser.add_argument('-n', '--increment_grid_search', help='Increment for gird search to find parameters with lowest'
                                                              'Euclidean distance. [0.02]', default=0.02,
                        type=float)
    parser.add_argument('-d', '--distance_cutoff', help='Max Euclidean distance to accept parameter combinations. '
                                                        '[0.01]', type=float, default=0.01)
    parser.add_argument('-g', '--n_generations', help='Number of generations since admixture. [15]', type=int,
                        default=15)
    parser.add_argument('-o', '--output', required=True, help='Output file name for figure showing box plots of '
                                                              'possible sex ratios for each ancestry.')
    parser.add_argument('--sf0', help='Initial female contribution. Result does not depend on this parameter '
                                      'if admixture happened more than 10 generations ago. [0.5]', type=float,
                        default=0.5)
    parser.add_argument('--sm0', help='Initial male contribution. Result does not depend on this parameter if admixture'
                                      'happened more than 10 generations ago. [0.5]', type=float, default=0.5)
    args = parser.parse_args()
    df, populations, ancestries = read_data(args.input)
    ancestries = sorted(ancestries)[::-1]
    # raw values
    fig, ax = plt.subplots(len(populations), 1, sharex=True, figsize=(10, 10))
    plt.subplots_adjust(hspace=0.65)
    base_colors = ['blue', 'orange', 'green'][::-1]
    colors = {anc: base_colors[i] if i < len(base_colors)
                                  else (np.random.random(), np.random.random(), np.random.random())
              for i, anc in enumerate(ancestries)}
    for n, population in enumerate(populations):
        for i, ancestry in enumerate(ancestries):
            parameter_estimates = constant_admixture(df.loc[population, '{}_A'.format(ancestry)],
                                                     df.loc[population, '{}_X'.format(ancestry)],
                                                     df.loc[population, 'pf'], args.increment_grid_search,
                                                     args.distance_cutoff, args.n_generations, args.sf0, args.sm0)
            sex_ratios = parameter_estimates[:, 0] / parameter_estimates[:, 1]
            sex_ratios = sex_ratios[(~np.isinf(sex_ratios)) & (~np.isnan(sex_ratios))]
            ax[n].boxplot(sex_ratios, positions=[i],
                          vert=False, widths=0.5, showfliers=False, whis=[2.5, 97.5],  medianprops={'color': 'black'})
            ax[n].scatter(sex_ratios, i + np.random.uniform(-0.15, 0.15, size=sex_ratios.shape[0]),
                       color=colors[ancestry], alpha=0.5, s=6, edgecolor='none')
        ax[n].set_yticks(range(0, len(ancestries)))
        ax[n].set_yticklabels(ancestries)
        ax[n].set_title(population)
    ax[n].set_xlabel(r"$s^f/s^m$")
    ax[n].set_xscale("log")
    ax[n].set_xlim([1e-2, 100])
    fig.savefig(args.output, dpi=500, bbox_inches='tight')


if __name__ == '__main__':
    main(sys.argv[1:])
