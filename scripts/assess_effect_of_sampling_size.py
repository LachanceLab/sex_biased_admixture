#!/usr/bin/env python
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from sklearn.utils import resample
from summarize_results import read_results, assign_k_to_ancestry, filter_individuals
from estimate_sex_ratios_single_pulse import calculate_sex_contributions


def main(argv):
    parser = argparse.ArgumentParser(description='Assess effect of sampling on inferred sex ratios '
                                                 'incl. confidence intervals')
    parser.add_argument('-qa', '--admixture_results_autosomes', nargs='+', help='.Q file from Admixture for autosomes')
    parser.add_argument('-qx', '--admixture_results_chrx', nargs='+', help='.Q file from Admixture for X chromosome')
    parser.add_argument('-f', '--fam_file', nargs='+', help='*.fam file generated by plink to get individual IDs')
    parser.add_argument('-t', '--target_population', help='Population for which to assess effect of sample size '
                                                          'on sex ratios. [amr]', default="amr")
    parser.add_argument('-b', '--bootstraps', default=None, help='Number of bootstrap resamples. '
                                                                 'Default is the number of samples and max. 10000',
                        type=int)
    parser.add_argument('--simulated_sex_ratios', help='Simulated sex ratios for African, European, and East Asian '
                                                       'contributions to American population. [2, 0.5, 1.25]',
                        default=[2, 0.5, 1.25])
    parser.add_argument('-p', '--plot_dir', help='Directory where to save plots [./plots/]', default='./plots/')
    parser.add_argument('-o', '--output', help='Name of output file to save figure to. [effect_of_sample_size.pdf]',
                        default='effect_of_sample_size.pdf', required=False)
    args = parser.parse_args()
    if args.plot_dir.endswith('/'):
        plot_dir = args.plot_dir
    else:
        plot_dir = args.plot_dir + '/'
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)
    sample_sizes = []
    generations = np.arange(2, 16, 1).tolist()
    # equilibrium model
    generations.append('inf')
    fig, ax = plt.subplots()
    base_colors = ['blue', 'orange', 'green']
    for n, (admixture_results_autosomes, admixture_results_chrx, fam_file) in \
            enumerate(zip(args.admixture_results_autosomes, args.admixture_results_chrx, args.fam_file)):
        admixture_autosomes, admixture_x = read_results(admixture_results_autosomes, admixture_results_chrx,
                                                        fam_file)
        admixture_autosomes = assign_k_to_ancestry(admixture_autosomes)
        admixture_x = assign_k_to_ancestry(admixture_x)
        ancestries = [col for col in admixture_autosomes.columns if col != "pop" and col != "sex"]
        if n == 0:
            colors = {anc: base_colors[i] if i < len(base_colors)
                      else (np.random.random(), np.random.random(), np.random.random())
                      for i, anc in enumerate(ancestries)}
        admixture_autosomes = admixture_autosomes[admixture_autosomes["pop"] == args.target_population]
        admixture_x = admixture_x[admixture_x["pop"] == args.target_population]
        sample_sizes.append(admixture_autosomes.shape[0])
        # filter individuals based on ancestry_proportion
        #individuals_to_remove = filter_individuals(admixture_autosomes)
        #admixture_autosomes.drop(individuals_to_remove, inplace=True)
        #admixture_x.drop(individuals_to_remove, inplace=True)
        # join
        ancestry_proportions = admixture_autosomes.join(admixture_x.drop(["sex"], axis=1), lsuffix='_A',
                                                        rsuffix='_X')
        sample_names = ancestry_proportions.index.values
        c_ancestry_proportions = ancestry_proportions.loc[sample_names, :]
        if args.bootstraps is None:
            n_resamples = max([min([c_ancestry_proportions.shape[0], 10000]), 500])
        else:
            n_resamples = args.bootstraps
        iterations_df = pd.DataFrame(index=np.repeat(np.arange(0, n_resamples), len(generations)))
        for i in range(n_resamples):
            resampled = resample(sample_names)
            resampled_autosomes = c_ancestry_proportions.loc[resampled, [f"{anc}_A" for anc in ancestries]]
            resampled_x = c_ancestry_proportions.loc[resampled, [f"{anc}_X" for anc in ancestries]]
            mean_sex = c_ancestry_proportions.loc[resampled, "sex"].mean()
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
                iterations_df.loc[i, '{}_sf/sm'.format(ancestry)] = [
                    x.round(2) if np.abs(x) < 10 else x.round(1) if np.abs(x) > 10
                                                                    and np.abs(
                        x) < 100 else int(x.round(0)) for x in females_to_male]
                iterations_df.loc[i, '{}_sm/sf'.format(ancestry)] = [
                    x.round(2) if np.abs(x) < 10 else x.round(1) if np.abs(x) > 10
                                                                    and np.abs(
                        x) < 100 else int(x.round(0)) for x in males_to_female]
        # summarize
        mean_vals = iterations_df.groupby("G").mean()
        # 95% CI
        ci_low = iterations_df.groupby('G').quantile(0.025)
        ci_up = iterations_df.groupby('G').quantile(0.975)
        ci_up.columns = [f"{col}_ci_up" for col in ci_up.columns]
        summary_df = ci_low.join(mean_vals, lsuffix="_ci_low", rsuffix="_mean").join(ci_up)
        # plot sex ratio with CI
        for m, ancestry in enumerate(ancestries):
            ax.errorbar(n + 0.1 * (m - 1), summary_df.loc[15, f"{ancestry}_sf/sm_mean"],
                        yerr=np.array([[summary_df.loc[15, f"{ancestry}_sf/sm_mean"] -
                                        summary_df.loc[15, f"{ancestry}_sf/sm_ci_low"],
                                        summary_df.loc[15, f"{ancestry}_sf/sm_ci_up"] -
                                        summary_df.loc[15, f"{ancestry}_sf/sm_mean"]]]).T,
                        color=colors[ancestry], marker='o', zorder=-1)
            ax.scatter(n + 0.1 * (m - 1), args.simulated_sex_ratios[m], marker='x', color='black')
    # Figure formatting
    labels = [r"$S_1$" if anc == 'afr' else r"$S_2$" if anc == 'eur' else r"$S_3$" for anc in ancestries]
    legend_handles = [Line2D([0], [0], color=colors[anc], label=l) for anc, l in zip(ancestries, labels)]
    legend_handles.append(Line2D([0], [0], color='black', label=r'Simulated $s^f/s^m$', ls='', marker='x'))
    ax.legend(handles=legend_handles, bbox_to_anchor=(0.5, -.14), ncol=4, loc='upper center')
    ax.set_xticks(np.arange(0, len(sample_sizes)))
    ax.set_xticklabels(sample_sizes)
    ax.set_ylabel(r'$s^f/s^m$')
    ax.set_xlabel('Sample size')
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.set_ylim([0.2, 3.5])
    fig.savefig(f"{plot_dir}{args.output}", dpi=500,  bbox_inches='tight')


if __name__ == '__main__':
    main(sys.argv[1:])
