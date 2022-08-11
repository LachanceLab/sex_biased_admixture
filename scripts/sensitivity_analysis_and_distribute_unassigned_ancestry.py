#!/usr/bin/env python
import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Polygon
from estimate_sex_ratios_single_pulse import read_data_summary_statistics, calculate_sex_contributions


def sensitivity_analysis(population, ancestry, data_file, figure_path):
    """
    Plot sensitivity analysis to show divergence (Figure 1)
    :param population: str, population which to use as an example for the sensitivity analysis.
                            Autosomal ancestry proportion for specified population and ancestry will be used.
    :param ancestry: str, ancestry which to use as an example for the sensitivity analysis.
                          Autosomal ancestry proportion for specified population and ancestry will be used.
    :param data_file: str, path to data file summary statistics, i.e., mean ancestry proportions
    :param figure_path: str, Path where to save figure of sensitivity analysis to

    """
    # load data file
    df, populations, ancestries = read_data_summary_statistics(data_file)
    # extract autosomal ancestry proportion for example
    autosomal_ancestry = df.loc[population, ancestry.upper() + '_A']
    # extract corresponding observed X chromosomal ancestry proportion
    observed_x_ancestry = df.loc[population, ancestry.upper() + '_X']
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
    ax[0].plot(x_ancestry, sex_ratio, lw=1, zorder=-1)
    # plot sex ratio based on observed ancestry proportions
    ax[0].scatter(observed_x_ancestry, sf_sm(observed_x_ancestry, autosomal_ancestry), marker='^', color='black', zorder=1)
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
    ax[1].add_patch(Polygon(patch_neg_too_much_x, closed=True, fill=True, color='darkgray', label='model failure',
                            zorder=-1))
    ax[1].add_patch(Polygon(patch_female, closed=True, fill=True,
                            color=(4.535548065502659e-06, 0.5400411873079152, 0.2349801931202921),
                            label='female-biased', zorder=-1))
    ax[1].add_patch(Polygon(patch_male, closed=True, fill=True,
                            color=(0.16305460270310312, 0.42335679778955587, 0.9999895609181783), label='male-biased',
                            zorder=-1))
    ax[1].add_patch(Polygon(patch_neg_too_little_x, closed=True, fill=True, color='darkgray', zorder=-1))
    ax[1].plot([0, 1], [0, 1], ls='--', color='black', lw=0.75)
    # plot observed ancestry proportions into space
    for pop in populations:
        for anc in ancestries:
            if df.loc[pop, anc.upper() + "_X"] >= 0.05 and df.loc[pop, anc.upper() + "_A"] >= 0.05:
                ax[1].scatter(df.loc[pop, anc.upper() + "_X"], df.loc[pop, anc.upper() + "_A"], marker='^',
                              color="black", s=6)
    ax[1].set_xlim([0, 1.0])
    ax[1].set_ylim([0, 1.0])
    ax[1].set_xlabel(r'X chromosomal ancestry proportion $\left(H^X_1\right)$', labelpad=10, fontdict=dict(fontsize=7))
    ax[1].set_ylabel(r'Autosomal ancestry proportion $\left(H^A_1\right)$', fontdict=dict(fontsize=7))
    ax[1].tick_params(labelsize=6)
    ax[1].set_aspect('equal')

    ax[1].legend(bbox_to_anchor=(1.025, 0.5), loc='center left', fontsize=7)
    ax[1].text(-0.22, ax[1].get_ylim()[1] + ax[1].get_ylim()[1] * 0.04, 'B', fontsize=12, weight='bold')

    fig.savefig(figure_path, dpi=500, bbox_inches='tight')
    plt.close()


def get_possible_ancestry_proportions(df, population, ancestries):
    """
    Get all possible ancestry proportions given estimated proportions and amount of unassigned ancestries
    :param df: pd.Dataframe, ancestry estimates
    :param population: str, current population of ancestries
    :param ancestries: list, ancestries
    :return: np.array, (n, 3) possible combinations of ancestry proportions
    """
    # calculate amount of missing ancestries
    missing_autosomal_ancestries = np.round(1 - np.round(df.loc[population,
                                                                [f'{ancestry}_A' for ancestry in ancestries]].values,
                                                         2).sum(), 2)
    missing_x_ancestries = np.round(1 - np.round(df.loc[population,
                                                        [f'{ancestry}_X' for ancestry in ancestries]].values,
                                                 2).sum(), 2)
    # create steps in which it can be distributed
    steps_autosomal = np.concatenate([np.arange(0, missing_autosomal_ancestries, 0.01),
                                      np.array([missing_autosomal_ancestries])])
    steps_x = np.concatenate([np.arange(0, missing_x_ancestries, 0.01), np.array([missing_x_ancestries])])
    # create all possible combination of distributing unassigned ancestries
    combinations_autosomal = np.stack(np.meshgrid(steps_autosomal, steps_autosomal, steps_autosomal), -1).reshape(-1, 3)
    combinations_x = np.stack(np.meshgrid(steps_x, steps_x, steps_x), -1).reshape(-1, 3)
    # select possible combination that distributed not less or more than is unassigned
    accepted_combinations_autosomal = combinations_autosomal[
        np.isclose(combinations_autosomal.sum(axis=1), missing_autosomal_ancestries)]
    accepted_combinations_x = combinations_x[np.isclose(combinations_x.sum(axis=1), missing_x_ancestries)]
    # calculate possible ancestry proportions
    possible_autosomal_proportions = accepted_combinations_autosomal + \
                                     np.round(df.loc[population, [f'{ancestry}_A' for ancestry in ancestries]].values,
                                              2)
    possible_x_proportions = accepted_combinations_x + \
                             np.round(df.loc[population, [f'{ancestry}_X' for ancestry in ancestries]].values, 2)
    # all possible combinations of autosomal and X chromosomal ancestry proportions
    parameter_space = np.stack([np.repeat(possible_autosomal_proportions, possible_x_proportions.shape[0], axis=0),
                                np.tile(possible_x_proportions.T, possible_autosomal_proportions.shape[0]).T])
    return parameter_space


def explore_parameter_space(parameter_space, ancestries):
    """
    Plot sex ratios of ancestry proportions if unassigned ancestries are distributed in all possible ways
    :param parameter_space: np.array, possible ancestry proportions
    :param ancestries: list, ancestries
    :return dict: conceivable sex ratios and initial sex_ratio per ancestry

    """
    sf_sm = lambda x, a: (3 * x - 2 * a) / (4 * a - 3 * x)
    possible_sex_ratios_ancestries = {}

    for i, ancestry in enumerate(ancestries):
        # get unique combinations of ancestry proportions
        proportion_combinations = pd.DataFrame.from_dict({'X': parameter_space[1, :, i],
                                                          'A': parameter_space[0, :, i]}).drop_duplicates()
        # sort ascending
        proportion_combinations.sort_values(['X', "A"], inplace=True)
        # compute sex ratios for all possible combinations
        x_proportions = proportion_combinations.loc[:, "X"]
        a_proportions = proportion_combinations.loc[:, "A"]
        possible_sex_ratios = sf_sm(x_proportions, a_proportions).values
        # flatten
        possible_sex_ratios = possible_sex_ratios.flatten()
        # set np.inf to nan
        possible_sex_ratios = np.where(possible_sex_ratios == np.inf, np.nan, possible_sex_ratios)
        # set misspecifications to nan, else take log
        possible_sex_ratios = np.where(possible_sex_ratios > 0, possible_sex_ratios, np.nan)
        possible_sex_ratios = possible_sex_ratios[~np.isnan(possible_sex_ratios)]

        possible_sex_ratios_ancestries[ancestry] = {'possible_values': possible_sex_ratios}

    return possible_sex_ratios_ancestries


def plot_range_of_possible_sex_ratios_target_populations(sex_ratios, populations_of_interest, ancestries, figure_path):
    """
    :param sex_ratios: dict, possible containing possible sex ratios after distributing unassigned ancestries in all
                             permissible ways and initially predicted sex ratio for each population of interest and ancestry
    :param populations_of_interest: list, populations for which to merge figures of exploration of parameter space
    :param ancestries: list, ancestries
    :param figure_path: str, Path where to save figure of box plots of possible sex ratios to
    """
    fig, ax = plt.subplots(len(populations_of_interest), 1, sharex=True, sharey=True)
    plt.subplots_adjust(hspace=0.65)
    ancestries = sorted(ancestries)[::-1]
    base_colors = ['blue', 'orange', 'green'][::-1]
    colors = {anc: base_colors[i] if i < len(base_colors)
              else (np.random.random(), np.random.random(), np.random.random())
              for i, anc in enumerate(ancestries)}
    for i, pop in enumerate(populations_of_interest):
        for n, ancestry in enumerate(ancestries):
            ax[i].boxplot(sex_ratios[pop][ancestry]['possible_values'], positions=[n],
                          vert=False, widths=0.5, showfliers=False,
                          whis=[2.5, 97.5], medianprops={'color': "black"})
            ax[i].scatter(sex_ratios[pop][ancestry]['possible_values'],
                          n + np.random.uniform(-0.15, 0.15, size=sex_ratios[pop][ancestry]['possible_values'].shape[0]),
                          color=colors[ancestry], alpha=0.5, s=6, edgecolor='none')
            # if not np.isnan(sex_ratios[pop][ancestry]['predicted_ratio']):
            ax[i].scatter(sex_ratios[pop][ancestry]['predicted_ratio'], n, marker='^', color='black', )

        ax[i].set_yticks(range(0, len(ancestries)))
        ax[i].set_yticklabels(ancestries)
        ax[i].set_title(pop)
    ax[i].set_xlabel(r"$s^f/s^m$")
    ax[i].set_xscale("log")
    ax[i].set_xlim([1e-2, 100])
    fig.savefig(figure_path, dpi=500, bbox_inches='tight')


def compute_sex_ratios_from_summary_statistics(input_file, output_file, figure_path, populations_of_interest=[]):
    """
    Compute sex ratios from summary statistics (mean ancestry proportions)
    :param input_file: str, Path to file with mean ancestry proportions
    :param output_file: str, Path where to write results to
    :param figure_path: str, Path where to save figure of box plots of possible sex ratios to
    :param populations_of_interest: list, populations for which to merge figures of exploration of parameter space
    """
    df, populations, ancestries = read_data_summary_statistics(input_file)

    generations = np.arange(2, 16, 1).tolist()
    generations.append('inf')
    complete_df = pd.DataFrame()
    possible_sex_ratios_pop_of_interest = dict()
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
            # handle infinities due to zero division --> replace with 999
            females_to_male = np.nan_to_num(females_to_male, posinf=999)
            males_to_female = np.nan_to_num(males_to_female, posinf=999)
            # round sex ratios to three significant digits
            c_df['{}_sf/sm'.format(ancestry)] = [x.round(2) if np.abs(x) < 10 else x.round(1) if np.abs(x) > 10
                                                                                                 and np.abs(
                x) < 100 else int(x.round(0)) for x in females_to_male]
            c_df['{}_sm/sf'.format(ancestry)] = [x.round(2) if np.abs(x) < 10 else x.round(1) if np.abs(x) > 10
                                                                                                 and np.abs(
                x) < 100 else int(x.round(0)) for x in males_to_female]
        complete_df = complete_df.append(c_df)
        if len(ancestries) != 3 or (np.isclose(df.loc[population, [f"{anc}_A" for anc in ancestries]].sum(), 1)
                                    and np.isclose(df.loc[population, [f"{anc}_X" for anc in ancestries]].sum(), 1)):
            continue
        parameter_space = get_possible_ancestry_proportions(df, population, ancestries)
        possible_sex_ratios_pop = explore_parameter_space(parameter_space,  ancestries)
        for ancestry in ancestries:
            # store initially predicted sex ratio
            possible_sex_ratios_pop[ancestry]['predicted_ratio'] = c_df['{}_sf/sm'.format(ancestry)].values[-1]

        if population in populations_of_interest:
            possible_sex_ratios_pop_of_interest[population] = possible_sex_ratios_pop

    # save to excel
    complete_df.reset_index(drop=True, inplace=True)
    complete_df.to_excel(output_file, index=False)
    if len(populations_of_interest) > 0:
        plot_range_of_possible_sex_ratios_target_populations(possible_sex_ratios_pop_of_interest,
                                                             populations_of_interest, ancestries,
                                                             figure_path)


def main(argv):
    parser = argparse.ArgumentParser(description="Re-implementations of Goldberg's recursion equation "
                                                 "to calculate the sex ratios for a single pulse admixture "
                                                 "https://www.genetics.org/content/201/1/263 "
                                                 "Provided data file contains ancestry estimates from"
                                                 "Micheletti et al. (2020)"
                                                 "https://doi.org/10.1016/j.ajhg.2020.06.012 "
                                                 "Sex ratios are "
                                                 "inferred based on the provided data first,"
                                                 "and then after distributing unassigned ancestries in all possible"
                                                 " ways such that the ancestries sum up to 100%. "
                                                 "Additionally, a sensitivity analysis of the models is conducted")
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
    parser.add_argument('-p1', '--populations_of_interest',
                        help="Populations for which to merge exploration of parameter space into one Figure", nargs='+',
                        required=False, default=[])
    parser.add_argument('-p2', '--population_of_interest_sensitivity_analysis',
                        help='Population which to use as an example for the sensitivity analysis. '
                             'It will use the autosomal ancestry proportion in the specified population '
                             'of the ancestry specified with -a  [Central America]',
                        required=False, default="Central America")
    parser.add_argument('-a', '--ancestry_of_interest_sensitivity_analysis',
                        help='Ancestry which to use as an example for the sensitivity analysis. '
                             'It will use the autosomal ancestry proportion in the population specified with -p2 '
                             'of the specified ancestry  [Afr]',
                        required=False, default="Afr")
    parser.add_argument('-o', '--output', help='Output filename. Will be an excel file.', required=True)
    parser.add_argument('--figure', required=True, help='Path to save figure with box plots of possible sex ratios to.')
    parser.add_argument('--figure_sensitivity', required=True, help='Path to save figure of senstivity analysis to.')
    args = parser.parse_args()
    sensitivity_analysis(args.population_of_interest_sensitivity_analysis,
                         args.ancestry_of_interest_sensitivity_analysis, args.input, args.figure_sensitivity)
    compute_sex_ratios_from_summary_statistics(args.input, args.output, args.figure, args.populations_of_interest)


if __name__ == '__main__':
    main(sys.argv[1:])
