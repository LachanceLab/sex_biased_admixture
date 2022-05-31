#!/usr/bin/env python
import sys
import argparse
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
from matplotlib.lines import Line2D
import logging
from sklearn.utils import resample


def read_results(admixture_results_autosomes, admixture_results_chrx, fam_file):
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
    # add sample populations
    populations = np.array([re.sub('[0-9]*', '', ind) for ind in indv_names_and_sex.iloc[:, 0]])
    admixture_autosomes.loc[:, 'pop'] = populations
    admixture_x.loc[:, 'pop'] = populations
    return admixture_autosomes, admixture_x


def assign_k_to_ancestry(df, ancestries=['afr', 'eur', 'ea']):
    """
    Identify K corresponding to source ancestries
    :param df: pd.DataFrame, admixture results
    :param ancestries: list, source populations
    :return: pd.DataFrame, with renamed columns according to main ancestries
    """
    mean_per_k = df.groupby("pop").mean().iloc[:, :-1]
    mean_per_k = mean_per_k.loc[ancestries]
    max_ind_k = np.argmax(mean_per_k.values, axis=0)
    if np.unique(max_ind_k).shape[0] < len(ancestries):
        logging.debug(f"Inadequate number of source populations inferred."
                      f"The mean contributions of source populations to different cluster were\n{mean_per_k}")
        raise ValueError('Not all source ancestries are represented. You likely used an inadequate number for K. '
                         'Check log file for more details')
    new_col_names = {col: ancestries[max_ind_k[i]] for i, col in enumerate(mean_per_k.columns)}
    df.rename(columns=new_col_names, inplace=True)
    new_df = pd.DataFrame(index=df.index)
    for anc in ancestries:
        try:
            new_df.loc[:, anc] = df.loc[:, anc].sum(axis=1)
        except ValueError:
            new_df.loc[:, anc] = df.loc[:, anc]
    new_df.loc[:, 'pop'] = df.loc[:, 'pop'].copy()
    new_df.loc[:, 'sex'] = df.loc[:, 'sex'].copy()
    return new_df


def plot_autosome_proportions_vs_x_proportions(proportions_autosomes, proportions_x, output_base):
    """
    Plot autosomal vs. X chromosomal ancestry proportions
    :param proportions_autosomes: pd.DataFrame, autosomal admixture results
    :param proportions_x: pd.DataFrame, X chromosomal admixture results
    :param output_base: str, output base to save plot per population to
    """
    np.random.seed(13)
    populations = proportions_autosomes.loc[:, 'pop'].unique()
    source_pop = np.sort(proportions_autosomes.columns[:-2])
    colors = {col: (np.random.random(), np.random.random(), np.random.random()) for col in source_pop}
    for pop in populations:
        fig, ax = plt.subplots()
        c_pop_autosomes = proportions_autosomes[proportions_autosomes['pop'] == pop]
        c_pop_x = proportions_x[proportions_x['pop'] == pop]
        for s_pop in source_pop:
            ax.scatter(c_pop_x[s_pop].values, c_pop_autosomes[s_pop].values, label=s_pop,
                       color=colors[s_pop], alpha=0.5)
        ax.set_xlabel('X chromosome proportion')
        ax.legend(bbox_to_anchor=(0.5, -.14), ncol=3, loc='upper center')
        ax.set_ylabel('Autosome proportion')
        ax.plot([0, 1], [0, 1], ls='--', c='black')
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        fig.savefig(output_base + f'_A_vs_X_{pop}.png', bbox_inches='tight')
        plt.close()


def plot_histogram_of_ancestry_proportions(proportions_autosomes, proportions_x, output_base):
    """
    Plot cumulative histograms of autosomal and X chromosomal ancestry proportions
    :param proportions_autosomes: pd.DataFrame, autosomal admixture results
    :param proportions_x: pd.DataFrame, X chromosomal admixture results
    :param output_base: str, output base to save plot per population to
    """
    np.random.seed(13)
    populations = proportions_autosomes.loc[:, 'pop'].unique()
    source_pop = np.sort(proportions_autosomes.columns[:-2])
    colors = {col: (np.random.random(), np.random.random(), np.random.random()) for col in source_pop}
    bins = np.arange(0, 1.01, 0.01)
    for pop in populations:
        fig, ax = plt.subplots()
        c_pop_autosomes = proportions_autosomes[proportions_autosomes['pop'] == pop]
        c_pop_x = proportions_x[proportions_x['pop'] == pop]
        for s_pop in source_pop:
            ax.hist(c_pop_x[s_pop].values, color=colors[s_pop], cumulative=True, ls='--', bins=bins,
                    histtype='step',
                    weights=np.repeat(1 / c_pop_x[s_pop].values.shape[0], c_pop_x[s_pop].values.shape[0]))
            ax.hist(c_pop_autosomes[s_pop].values, color=colors[s_pop], cumulative=True, ls='-', bins=bins,
                    histtype='step',
                    weights=np.repeat(1 / c_pop_autosomes[s_pop].values.shape[0],
                                      c_pop_autosomes[s_pop].values.shape[0]))
        ax.set_xlabel('Ancestry proportion')
        ax.set_ylabel('Cumulative density')
        handles = [Line2D([0], [0], color=colors[source_pop[n]], ls='--') for n in range(len(source_pop))]
        handles.extend([Line2D([0], [0], color='black', ls='-'), Line2D([0], [0], color='black', ls='--')])
        labels = source_pop.copy().tolist()
        labels.extend(['Autosomes', 'X chromosomes'])
        ax.legend(handles, labels, bbox_to_anchor=(0.5, -.14), loc='upper center', ncol=3)
        fig.savefig(output_base + f"_hist_a_x_proportions_{pop}.png", bbox_inches='tight')


def plot_admixture_results(df, outfile):
    """
    Plot admixture results
    :param df: pd.DataFrame, admixture results
    :param output_base: str, output base to save plot per population to
    """
    np.random.seed(13)
    populations = df['pop'].unique()
    source_pop = np.sort(df.columns[:-2])
    colors = {col: (np.random.random(), np.random.random(), np.random.random()) for col in source_pop}
    fig, ax = plt.subplots()
    proportions_dict = {}
    colors_dict = {}
    for i in range(len(source_pop)):
        proportions_dict[i] = np.zeros(df.shape[0])
        colors_dict[i] = np.zeros((df.shape[0], 3))
    i = 0
    n_samples = []
    for pop in populations:
        c_df = df[df['pop'] == pop]
        proportions = c_df.loc[:, source_pop]
        median_prop = np.median(proportions, axis=0)
        priorities = source_pop[np.argsort(median_prop)][::-1].tolist()
        for n in range(len(source_pop)):
            proportions_dict[n][i: i + c_df.shape[0]] = c_df.loc[:, priorities[n]]
            colors_dict[n][i: i + c_df.shape[0]] = colors[priorities[n]]
        i += c_df.shape[0]
        n_samples.append(c_df.shape[0])
    bottom = np.zeros(df.shape[0])
    for n in range(len(source_pop)):
        ax.bar(np.arange(0, df.shape[0]), proportions_dict[n], color=colors_dict[n], bottom=bottom, width=1)
        bottom += proportions_dict[n]
    n_samples = np.array(n_samples)
    cum_samples = np.cumsum(n_samples)
    n_samples = n_samples / 2
    xticks = cum_samples - n_samples
    ax.set_xticks(xticks)
    ax.set_ylim([0, 1.0])
    ax.set_xlim([0, df.shape[0]])
    ax.set_xticklabels([pop.upper() for pop in populations])
    ax.set_ylabel('Admixture proportion')
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    fig.savefig(outfile, bbox_inches='tight')
    plt.close()


def write_mean_ancestry_proportions(df, output_file):
    """
    Write mean ancestry proportions to file
    :param df: pd.DataFrame, mean ancestry proportions of source populations for each sample populations
    :param output_file: str, Path to file where to write proportions to
    """
    # bringing columns into the right order
    columns = {col: col.replace('ea', 'na').replace('_A', '_ZA').upper() for col in df.columns}
    df.rename(columns=columns, inplace=True)
    columns = sorted([col for col in df.columns])
    df = df.loc[:, columns]
    columns = {col: col.replace('_ZA', '_A').replace('SEX', 'pf') for col in df.columns}
    df.rename(columns=columns, inplace=True)
    df.index.name = ''
    df.to_csv(output_file, header=True, index=True)


def filter_individuals(admixture_autosomes):
    """
    Filtered out African/European/East Asian individuals with less than 95% African/Eurpean/East Asian ancestry as well
    as admixed individuals with more than 95% of one ancestry
    :param admixture_autosomes: pd.DataFrame, Admixture results
    :return: np.array, IIDs of individuals that not passed filtering
    """
    # filter Africans with less 95% African ancestry
    afr_to_remove = admixture_autosomes[(admixture_autosomes['pop'] == 'afr') &
                                        (admixture_autosomes.afr < 0.95)].index.values
    logging.info('Removing {} African individuals that have <95% African ancestry'.format(afr_to_remove.shape[0]))
    # filter Europeans with less 95% European ancestry
    eur_to_remove = admixture_autosomes[(admixture_autosomes['pop'] == 'eur') &
                                        (admixture_autosomes.eur < 0.95)].index.values
    logging.info('Removing {} European individuals that have <95% European ancestry'.format(eur_to_remove.shape[0]))
    # filter East Asians with less 95% East Asian ancestry
    ea_to_remove = admixture_autosomes[(admixture_autosomes['pop'] == 'ea') &
                                       (admixture_autosomes.ea < 0.95)].index.values
    logging.info('Removing {} East Asian individuals that have <95% East Asian ancestry'.format(ea_to_remove.shape[0]))
    # filter Americans with less 5% African ancestry
    amr_to_remove = admixture_autosomes[(admixture_autosomes['pop'] == 'amr') &
                                        (admixture_autosomes.loc[:, ['afr', 'eur', 'ea']].max(axis=1) > 0.95)].index.values
    logging.info('Removing {} American individuals that have >95% of one ancestry'.format(amr_to_remove.shape[0]))
    individuals_to_remove = np.concatenate([amr_to_remove, eur_to_remove, ea_to_remove, afr_to_remove])
    return individuals_to_remove


def compute_haplogroup_frequencies(df):
    """
    Compute haplogroup frequencies
    :param df: pd.DataFrame, haplogroup assignments
    :return: pd.DataFrame, haplogroup frequencies of sampled populations
    """
    haplogroups = np.sort(df.haplogroup.unique())
    populations = df.population.unique()
    frequencies = pd.DataFrame(index=populations, columns=haplogroups)
    for pop in populations:
        freqs = np.zeros_like(haplogroups)
        c_df = df[df.population == pop]
        counts = c_df.haplogroup.value_counts()
        for hg, c in zip(counts.index, counts):
            freqs[haplogroups == hg] = c
        freqs = freqs / freqs.sum()
        frequencies.loc[pop, :] = freqs
    return frequencies


def bootstrap_results(ancestry_proportions_haplogroups, n_resamples=10000):
    """
    Bootstrap individuals to compute confidence intervals for summary statistics (mean ancestry proportions,
    and sex ratios)
    :param ancestry_proportions_haplogroups: pd.DataFrame, individual level autosomal and X chromosomal ancestry
                                                           proportions as well as mtDNA and Y chromosome haplogroups
    :param n_resamples: int, number of resampling steps to perform
    :return: pd.DataFrame, Mean summary statistics with confidence intervals
    """
    ancestries = [col.split('_')[0] for col in ancestry_proportions_haplogroups.columns if col.endswith('A')]
    populations = ancestry_proportions_haplogroups.loc[:, "pop"].unique()
    sf_sm = lambda x, a: (3 * x - 2 * a) / (4 * a - 3 * x)

    bootstrapped_results = []
    for pop in populations:
        c_results = ancestry_proportions_haplogroups.loc[ancestry_proportions_haplogroups.loc[:, "pop"] == pop]
        sample_names = c_results.index.values
        if n_resamples is None:
            c_resamples = min([c_results.shape[0], 10000])
        iterations_df = pd.DataFrame(index=np.arange(0, c_resamples))
        for i in range(c_resamples):
            resampled = resample(sample_names)
            y_haplogroups = c_results.loc[resampled, "haplogroup_Y"].copy().dropna()
            mt_haplogroups = c_results.loc[resampled, 'haplogroup_mt']
            resampled_autosomes = c_results.loc[resampled, [f"{anc}_A" for anc in ancestries]]
            resampled_x = c_results.loc[resampled, [f"{anc}_X" for anc in ancestries]]
            for anc in ancestries:
                mt_haplo_freq = mt_haplogroups[mt_haplogroups == anc].shape[0] / mt_haplogroups.shape[0]
                y_haplo_freq = y_haplogroups[y_haplogroups == anc].shape[0] / y_haplogroups.shape[0]
                try:
                    haplogroup_imbalance = mt_haplo_freq / y_haplo_freq
                except ZeroDivisionError:
                    haplogroup_imbalance = np.inf
                mean_autosomal_prop = resampled_autosomes.loc[:, f"{anc}_A"].mean()
                mean_x_prop = resampled_x.loc[:, f"{anc}_X"].mean()
                sex_ratio = sf_sm(mean_x_prop, mean_autosomal_prop)
                iterations_df.loc[i, f'{anc}_A'] = mean_autosomal_prop
                iterations_df.loc[i, f'{anc}_X'] = mean_x_prop
                iterations_df.loc[i, f'{anc}_sf_sm'] = sex_ratio
                iterations_df.loc[i, f'{anc}_mt'] = mt_haplo_freq
                iterations_df.loc[i, f'{anc}_Y'] = y_haplo_freq
                iterations_df.loc[i, f'{anc}_mt_Y'] = haplogroup_imbalance

        mean_vals = iterations_df.mean(axis=0)
        mean_vals.index = [f"mean_{ind}" for ind in mean_vals.index.values]
        ci_low = iterations_df.quantile(0.025, axis=0)
        ci_low.index = [f"ci_low_{ind}" for ind in ci_low.index.values]
        ci_up = iterations_df.quantile(0.975, axis=0)
        ci_up.index = [f"ci_up_{ind}" for ind in ci_up.index.values]

        summary_df = pd.concat([ci_low, mean_vals, ci_up]).to_frame(pop)
        bootstrapped_results.append(summary_df.T)
    bootstrapped_results = pd.concat(bootstrapped_results)
    return bootstrapped_results


def main(argv):
    parser = argparse.ArgumentParser(description="Summarize simulation results, i.e., sex ratios, "
                                                 "admixture proportions, haplogroup frequencies, "
                                                 "and haplogroup imbalances including confidence intervals.")
    parser.add_argument('-qa', '--admixture_results_autosomes', nargs='+', help='.Q file from Admixture for autosomes')
    parser.add_argument('-qx', '--admixture_results_chrx', nargs='+', help='.Q file from Admixture for X chromosome')
    parser.add_argument('-hm', '--haplogroups_mtDNA', nargs='+', help='Haplogroup assignments mtDNA')
    parser.add_argument('-hy', '--haplogroups_Y', nargs='+', help='Haplogroup assignments Y chromosomes')
    parser.add_argument('-f', '--fam_file', nargs='+', help='*.fam file generated by plink to get individual IDs')
    parser.add_argument('-p', '--plot_dir', help='Directory where to save plots [./plots/]', default='./plots/')
    parser.add_argument('-rax', '--results_a_x', nargs='+', help='File to write autosomal and X chromosomal '
                                                                 'ancestry proportions to')
    parser.add_argument('-rmy', '--results_mt_y', nargs='+', help='File to write mtDNA and Y chromosomal haplogroup '
                                                                  'frequencies to')
    parser.add_argument('-b', '--bootstraps', default=None, help='Number of bootstrap resamples. Default is the '
                                                                 'number of samples and max. 10000', type=int)
    parser.add_argument('-rb', '--results_bootstrapped', nargs="+", help='File to write bootstrapped summary '
                                                                         'statistics to')
    parser.add_argument('--log', help='Log file [results/summarizing_results.log]',
                        default='results/summarizing_results.log')
    args = parser.parse_args()
    logging.basicConfig(filename=args.log, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S')
    logger = logging.getLogger()
    if args.plot_dir.endswith('/'):
        plot_dir = args.plot_dir
    else:
        plot_dir = args.plot_dir + '/'
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)
    results_dir_hg_freq = '/'.join(args.results_mt_y[0].split('/')[:-1])
    if not os.path.isdir(results_dir_hg_freq):
        os.makedirs(results_dir_hg_freq)
    results_dir_ap_freq = '/'.join(args.results_a_x[0].split('/')[:-1])
    if not os.path.isdir(results_dir_ap_freq):
        os.makedirs(results_dir_ap_freq)
    for admixture_results_autosomes, admixture_results_chrx, fam_file, \
        haplogroups_mtDNA, haplogroups_Y, results_a_x,\
        results_mt_y, results_bootstrapped in zip(args.admixture_results_autosomes,
                                                  args.admixture_results_chrx,
                                                  args.fam_file,
                                                  args.haplogroups_mtDNA,
                                                  args.haplogroups_Y,
                                                  args.results_a_x,
                                                  args.results_mt_y,
                                                  args.results_bootstrapped):
        logging.info('Parsing results {}'.format(admixture_results_autosomes.split('_chr1.3.Q')[0]))
        admixture_autosomes, admixture_x = read_results(admixture_results_autosomes, admixture_results_chrx, fam_file)
        haplogroups_mtDNA = pd.read_csv(haplogroups_mtDNA, sep='\t', index_col=0)
        haplogroups_Y = pd.read_csv(haplogroups_Y, sep='\t', index_col=0)
        populations_mtDNA = np.array([re.sub('[0-9]*', '', ind) for ind in haplogroups_mtDNA.index.values])
        populations_Y = np.array([re.sub('[0-9]*', '', ind) for ind in haplogroups_Y.index.values])
        haplogroups_mtDNA.loc[:, 'population'] = populations_mtDNA
        haplogroups_Y.loc[:, 'population'] = populations_Y
        admixture_autosomes = assign_k_to_ancestry(admixture_autosomes)
        admixture_x = assign_k_to_ancestry(admixture_x)
        # remove individuals that did not pass filtering
        logging.info('Filtering individuals')
        # filter individuals based on average ancestry_proportion 
        # assumes equal length of autosomes and X chromosome
        #average_ancestry_proportions = admixture_autosomes.copy()
        #average_ancestry_proportions.loc[:, ["afr", "eur", "ea"]] += admixture_x.loc[:, ["afr", "eur", "ea"]]
        #average_ancestry_proportions.loc[:, ["afr", "eur", "ea"]] /= 2
        individuals_to_remove = filter_individuals(admixture_autosomes)
        admixture_autosomes.drop(individuals_to_remove, inplace=True)
        admixture_x.drop(individuals_to_remove, inplace=True)
        haplogroups_mtDNA.drop(individuals_to_remove, inplace=True)
        haplogroups_Y.drop(individuals_to_remove[np.isin(individuals_to_remove, haplogroups_Y.index)], inplace=True)

        # visualize ADMIXTURE results
        logging.info('Visualizing ADMIXTURE results')
        plot_admixture_results(admixture_autosomes,
                               f'{plot_dir}{admixture_results_autosomes.split("/")[-1].replace(".3.Q", "_admixture.png")}')
        plot_admixture_results(admixture_x,
                               f'{plot_dir}{admixture_results_chrx.split("/")[-1].replace(".3.Q", "_admixture.png")}')
        plot_autosome_proportions_vs_x_proportions(admixture_autosomes, admixture_x,
                                                   f'{plot_dir}{admixture_results_chrx.split("/")[-1].replace(".3.Q", "")}')
        plot_histogram_of_ancestry_proportions(admixture_autosomes, admixture_x,
                                               f'{plot_dir}{admixture_results_chrx.split("/")[-1].replace(".3.Q", "")}')
        # join results
        ancestry_proportions = admixture_autosomes.join(admixture_x.drop(["pop", "sex"], axis=1),
                                                        lsuffix='_A', rsuffix='_X')
        logging.info('Computing mean ancestry proportions')
        # take average per populations
        mean_ancestry_proportions = ancestry_proportions.groupby('pop').mean()
        logging.info(f'Writing mean ancestry proportions to {results_a_x}')
        # write results
        write_mean_ancestry_proportions(mean_ancestry_proportions, results_a_x)
        logging.info("Computing mean haplogroup frequencies")
        # compute haplogroup frequencies
        freq_mtDNA = compute_haplogroup_frequencies(haplogroups_mtDNA)
        freq_Y = compute_haplogroup_frequencies(haplogroups_Y)
        # merge
        haplogroup_frequencies = freq_Y.join(freq_mtDNA, lsuffix='_Y', rsuffix='_MT')
        columns = {col: col.replace("ea", 'na').upper() for col in haplogroup_frequencies.columns}
        haplogroup_frequencies.rename(columns=columns, inplace=True)
        logging.info(f"Writing mean haplogroup frequencies to {results_mt_y}")
        # save
        haplogroup_frequencies.to_csv(results_mt_y, header=True, index=True)
        # bootstrap results
        logging.info("Bootstrapping results")
        haplogroups = haplogroups_mtDNA.drop("population", axis=1).join(haplogroups_Y.drop("population", axis=1),
                                                                        lsuffix='_mt', rsuffix='_Y')
        ancestry_prop_and_haplogroups = ancestry_proportions.join(haplogroups)
        bootstrapped_results = bootstrap_results(ancestry_prop_and_haplogroups, args.bootstraps)
        logging.info(f"Writing bootstrap results to {results_bootstrapped}")
        bootstrapped_results.to_csv(results_bootstrapped, sep='\t', header=True, index=True)


if __name__ == '__main__':
    main(sys.argv[1:])
