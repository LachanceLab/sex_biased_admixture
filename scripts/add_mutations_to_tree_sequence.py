#!/usr/bin/env python
import pickle
import pyslim
import msprime
import sys
import argparse
import numpy as np
import os
import multiprocessing as mp
import matplotlib.pyplot as plt
import logging
import warnings
import re
from multiprocessing_logging import install_mp_handler


def tree_heights(ts):
    """
    Taken from pyslim manual
    :param ts: treeSequence
    :return: array, tree heights
    """
    # Calculate tree heights, giving uncoalesced sites the maximum time
    heights = np.zeros(ts.num_trees + 1)
    for tree in ts.trees():
        if tree.num_roots > 1:  # not fully coalesced
            heights[tree.index] = ts.metadata['SLiM']['generation']
        else:
            children = tree.children(tree.root)
            real_root = tree.root if len(children) > 1 else children[0]
            heights[tree.index] = tree.time(real_root)
    heights[-1] = heights[-2]  # repeat the last entry for plotting with step
    return heights


def visualize_recapitation(pre_recap, post_recap):
    """
    Visualize effect of recapitation
    :param pre_recap: TreeSequence, ts before recapitation
    :param post_recap: TreeSequence, ts after recapitation
    """
    fig, ax = plt.subplots(2, 1, sharex=True, sharey=True)
    # Plot tree heights before recapitation
    breakpoints = list(pre_recap.breakpoints())
    heights = tree_heights(pre_recap)
    ax[0].step(breakpoints, heights, where='post')
    ax[0].set_title('Prior recapitation')
    ax[0].set_ylabel('TMRCA')
    # ax[0].set_xlabel('Position')
    max_prior = heights.max()
    # Plot the tree heights after recapitation
    breakpoints = list(post_recap.breakpoints())
    heights = tree_heights(post_recap)
    ax[1].step(breakpoints, heights, where='post')
    ax[1].axhline(max_prior, color='black', ls='--')
    ax[1].set_xlabel('Position')
    ax[1].set_ylabel('TMRCA')
    ax[1].set_title('Post recapitation')
    fig.savefig('plots/vis_recapitation.png', bbox_inches='tight')
    plt.close()


def get_sfs(ts, pop, pop_index, filename, chr):
    """
    Plot site frequency spectrum
    :param ts: TreeSequence
    :param pop: str, current population
    :param filename: str, path to save figure to
    :param chr: str, current chromosome
    """
    fig, ax = plt.subplots()
    bins = np.arange(0, 101, 1)
    # plot real data
    if pop in ['afr', 'eur', 'ea']:
        if pop == 'afr':
            pop = 'yri'
        elif pop == 'eur':
            pop = 'ceu'
        elif pop == 'ea':
            pop = 'chb'
        if re.match('chr[0-9]+', chr):
            real_data = np.loadtxt(f'data/daf_frequencies_kg/allele_frequencies_{pop}_autosomes.txt')
        elif chr == 'chrX':
            real_data = np.loadtxt(f'data/daf_frequencies_kg/allele_frequencies_{pop}_chrX.txt')
        elif chr == 'chrY':
            real_data = np.loadtxt(f'data/daf_frequencies_kg/allele_frequencies_{pop}_chrY.txt')
        elif chr == 'mtDNA':
            real_data = np.loadtxt(f'data/daf_frequencies_kg/allele_frequencies_{pop}_mtDNA.txt')
        real_data = real_data[real_data > 0]
        real_data *= 100
        hist, bin_edges = np.histogram(real_data, bins=bins, weights=np.repeat(1 / real_data.shape[0],
                                                                               real_data.shape[0]))
        ax.scatter(bin_edges[:-1], hist, marker='^', color='blue', label='real')
    corresponding_nodes = ts.samples(pop_index)
    ts = ts.simplify(samples=corresponding_nodes)
    afs = ts.allele_frequency_spectrum(span_normalise=False, polarised=True)
    afs = afs[1:]
    frequencies = np.ones_like(afs)
    frequencies = np.cumsum(frequencies) / frequencies.shape[0] * 100
    afs /= afs.sum()
    hist, bin_edges = np.histogram(frequencies, bins=bins, weights=afs)
    ax.scatter(bin_edges[:-1], hist, marker='+', color='orange', label='simulations')

    ax.legend()
    fig.savefig(f"plots/{filename}", bbox_inches='tight')
    plt.close()


def add_mutations(tree_path, out_tree_path, Ne, mu, mtDNA_mu, chr_length, mtDNA_length, populations, recapitate, force):
    """
    Recapitate tree sequence, simply and overlay with mutations
    :param tree_path: str, path to input tree sequence
    :param out_tree_path: str, path were to dump mutated tree sequence
    :param Ne: float, ancestral population size
    :param mu: float, mutation rate
    :param mtDNA_mu: float, mutation rate mitochondrial DNA
    :param chr_length: int, chromosome length
    :param mtDNA_length: int, mtDNA length
    :param populations: list, acronyms of populations that are considered
    :param recapitate: boolean, whether to recaptitate or not
    :param force: boolean, whether to run full pipeline if out_tree_path already exists
    :return: TreeSequence, dict, tree sequence with mutations, pointer to nodes with Y-chromsomes and mtDNAs
    """
    # load tree sequence
    logging.info(f'Loading tree sequence {tree_path}')

    if os.path.isfile(out_tree_path) and not force:
        mutated = pyslim.load(out_tree_path)
        f = open(out_tree_path.replace('.trees', '.pkl'), 'rb')
        sample_genotype_mapping_by_population = pickle.load(f)
        f.close()

    else:
        ts = pyslim.load(tree_path)
        if recapitate:
            # recapitate (Ne = ancestral Ne)
            logging.info(f'Recapitating {tree_path}')
            ts_recap = pyslim.recapitate(ts, Ne)
            visualize_recapitation(ts, ts_recap)
            ts = ts_recap

        # simplify
        logging.info(f'Simplifying {tree_path}')
        ts.simplify()

        # assert all sites coalesced
        for t in ts.trees():
            if not t.num_roots == 1:
                warnings.warn("Tree {} did not completely coalesced."
                              "Consider setting --recapitate flag.".format(tree_path, t.interval[0], t.interval[1]))
                break
        sample_genotype_mapping_by_population = {}
        for i, pop in enumerate(populations):
            sample_genotype_mapping_by_population[pop] = {}
            node_ids = ts.samples(i + 1)
            sample_genotype_mapping_by_population[pop]['nodes'] = node_ids
            pop_var = [var for var in ts.variants(samples=node_ids)]
            Y_pop_var = [var for var in pop_var if var.site.position == 2 * chr_length - 1]
            mtDNA_var = [var for var in pop_var if var.site.position == 2 * chr_length + mtDNA_length - 1]
            # nodes that have Y chromosome
            Y_chroms = node_ids[np.where(Y_pop_var[0].genotypes == 1)[0]]
            sample_genotype_mapping_by_population[pop]['Y'] = Y_chroms
            # nodes that have maternal mtDNA
            mtDNA = node_ids[np.where(mtDNA_var[0].genotypes == 1)[0]]
            sample_genotype_mapping_by_population[pop]['mtDNA'] = mtDNA

        f = open(out_tree_path.replace('.trees', '.pkl'), 'wb')
        pickle.dump(sample_genotype_mapping_by_population, f)
        f.close()

        # add mutations
        logging.info(f'Overlaying mutations {tree_path}')
        mutation_rate_map = msprime.RateMap(position=[0, 2 * chr_length, 2 * chr_length + mtDNA_length],
                                            rate=[mu, mtDNA_mu])
        mutated = msprime.sim_mutations(ts, rate=mutation_rate_map, keep=False)
        mutated.dump(out_tree_path)
    return mutated, sample_genotype_mapping_by_population


def reverse_map_individuals(old_new_mapping, individuals_nodes_mapping, ts):
    """
    Map node IDs after simplifying to previous ones
    :param old_new_mapping: np.array, Has the shape (#nodes) where the uth entry gives the new node ID of
                                    uth node in the old mapping
    :param individuals_nodes_mapping: dict, maps old node ids to individual ids
    :param ts: TreeSequence, simplified tree sequence
    :return: list, list, old individual ids, new individual ids
    """
    old_individuals = []
    new_individuals = []
    for indv, nds in individuals_nodes_mapping.items():
        for nd in nds:
            new_node_id = old_new_mapping[nd]
            if new_node_id == -1:
                continue
            new_indv = ts.node(new_node_id).individual
            new_individuals.append(new_indv)
            old_individuals.append(indv)
            break
    return old_individuals, new_individuals 


def extract_population_specific_chromosomes(ts, pointers_Y_mtDNA, output_base, chr_length,
                                            mtDNA_length, populations, sample_size, target_population,
                                            off_target_sample_size, plot_sfs, threads):
    """
    Write different chromosomes to VCF for each population
    :param ts: TreeSequence, tree sequence to sample from
    :param pointers_Y_mtDNA: dict, pointers to nodes with mtDNA and Y chromosomes
    :param output_base: str, string base for output file
    :param chr_length: int, chromosome length
    :param mtDNA_length: int, length of mtDNA
    :param populations: list, population acronyms
    :param sample_size: int, number of individuals to sample from each population
    :param target_population: str, population for which to extract full sample size
    :param off_target_sample_size: int, sample size for off-target populations
    :param plot_sfs: boolean, whether to plot SFS or not
    :param threads: int, number of processes
    """
    # get current populations
    sampled_individuals_all_pop = []
    sampled_nodes_all_pop = []
    sampled_Y_all_pop = []
    sampled_X_all_pop = []
    sampled_mtDNA_all_pop = []
    for i, pop in enumerate(populations):
        logging.info(f'Extracting chromosome specific tree sequences for {pop}')
        # population specific nodes
        nodes = pointers_Y_mtDNA[pop]['nodes']
        # nodes with Y chromosomes
        Y_chromosomes = pointers_Y_mtDNA[pop]['Y']
        # nodes with X chromosomes
        X_chromosomes = nodes[~np.isin(nodes, Y_chromosomes)]
        # nodes with maternal mtDNA
        mtDNA = pointers_Y_mtDNA[pop]['mtDNA']
        # population specific individuals
        individuals = list(set([nd.individual for nd in ts.tables.nodes[nodes]]))
        if pop == target_population and len(individuals) < sample_size:
            logging.info(f'Population {pop} only has {len(individuals)} individuals. '
                         f'Reducing sample size for {pop} from {sample_size} to {len(individuals)}')
            c_sample_size = len(individuals)
        elif pop == target_population and len(individuals) >= sample_size:
            c_sample_size = sample_size
        elif pop != target_population and len(individuals) >= off_target_sample_size:
            c_sample_size = off_target_sample_size
        else:
            c_sample_size = len(individuals)
        # taken random sample of individuals
        sampled_individuals = np.random.choice(individuals, size=c_sample_size, replace=False)
        sampled_nodes = np.concatenate([ts.individual(indv).nodes for indv in sampled_individuals])
        sampled_Y = Y_chromosomes[np.isin(Y_chromosomes, sampled_nodes)]
        sampled_X = X_chromosomes[np.isin(X_chromosomes, sampled_nodes)]
        sampled_mtDNA = mtDNA[np.isin(mtDNA, sampled_nodes)]
        sampled_individuals_all_pop.append(sampled_individuals)
        sampled_nodes_all_pop.append(sampled_nodes)
        sampled_Y_all_pop.append(sampled_Y)
        sampled_X_all_pop.append(sampled_X)
        sampled_mtDNA_all_pop.append(sampled_mtDNA)

    sampled_indv_amr = sampled_individuals.copy()
    # merge samples from all populations
    sampled_individuals = np.concatenate(sampled_individuals_all_pop)
    sampled_nodes = np.concatenate(sampled_nodes_all_pop)
    sampled_Y = np.concatenate(sampled_Y_all_pop)
    sampled_X = np.concatenate(sampled_X_all_pop)
    sampled_mtDNA = np.concatenate(sampled_mtDNA_all_pop)
    sampled_individuals_node_mapping = {indv: ts.individual(indv).nodes for indv in sampled_individuals}

    # extract tree sequences corresponsing to chromosome 1
    ts_chr1, node_map_chr1 = ts.keep_intervals(np.array([[0, chr_length]]),
                                               simplify=False).simplify(sampled_nodes, map_nodes=True)
    # extract tree sequences corresponsing to sex chromosomes
    ts_Y_chromosomes, node_map_Y = ts.keep_intervals(np.array([[chr_length, 2 * chr_length]]),
                                                     simplify=False).simplify(sampled_Y, map_nodes=True)
    ts_X_chromosomes, node_map_X = ts.keep_intervals(np.array([[chr_length, 2 * chr_length]]),
                                                     simplify=False).simplify(sampled_X, map_nodes=True)
    # extract tree sequences for mtDNA
    ts_mtDNAs, node_map_mtDNA = ts.keep_intervals(np.array([[2 * chr_length, 2 * chr_length + mtDNA_length]]),
                                                  simplify=False).simplify(sampled_mtDNA, map_nodes=True)

    sampled_individuals, individuals_chr1 = reverse_map_individuals(node_map_chr1, sampled_individuals_node_mapping,
                                                                    ts_chr1)
    sampled_individuals_Y, individuals_Y = reverse_map_individuals(node_map_Y, sampled_individuals_node_mapping,
                                                                   ts_Y_chromosomes)
    sampled_individuals_X, individuals_X = reverse_map_individuals(node_map_X, sampled_individuals_node_mapping,
                                                                   ts_X_chromosomes)
    sampled_individuals_mtDNA, individuals_mtDNA = reverse_map_individuals(node_map_mtDNA,
                                                                           sampled_individuals_node_mapping, ts_mtDNAs)


    # get sex of individuals and write to file
    sampled_individuals_sex = ['F' if ts_chr1.individual(indv).metadata['sex'] == 0 else 'M'
                               for indv in individuals_chr1]
    sampled_individuals_population = [populations[ts_chr1.node(ts_chr1.individual(indv).nodes[0]).population]
                                      for indv in individuals_chr1]
    with open(f"{output_base}_sex.tab", 'w') as sex_pointer:
        sex_pointer.write('#IID\tsex\n')
        for indv, sex, pop in zip(sampled_individuals, sampled_individuals_sex, sampled_individuals_population):
            sex_pointer.write(f'{pop}{indv + 1}\t{sex}\n')
    sex_pointer.close()
    ## plot SFS
    if plot_sfs:
        for pop_index, pop in enumerate(populations):
            logging.info(f'Plotting SFS {output_base} for {pop}')
            for tree_seqs, chr in zip([ts_chr1, ts_X_chromosomes, ts_Y_chromosomes, ts_mtDNAs],
                                      ['chr1', 'chrX', 'chrY', 'mtDNA']):
                get_sfs(tree_seqs, pop, pop_index, f"{output_base.split('/')[-1]}_{pop}_{chr}_sfs.png", chr)

    # write .pop file for supervised admixture runs
    for ts, chr, individuals in zip([ts_chr1, ts_X_chromosomes, ts_Y_chromosomes, ts_mtDNAs],
                                    ['chr1', 'chrX', 'chrY', 'mtDNA'],
                                    [individuals_chr1, individuals_X, individuals_Y, individuals_mtDNA]):
        write_pop_file(ts, individuals, populations, target_population, f"{output_base}_{chr}.pop")

    ## write VCFs in parallel
    ready_to_map = [(tree_seqs, indvs, indv_names, f"{output_base}_{chr}.vcf", chr, populations) for
                    tree_seqs, indvs, indv_names, chr in zip([ts_chr1, ts_X_chromosomes, ts_Y_chromosomes, ts_mtDNAs],
                                                             [individuals_chr1, individuals_X, individuals_Y,
                                                              individuals_mtDNA],
                                                             [sampled_individuals, sampled_individuals_X,
                                                              sampled_individuals_Y, sampled_individuals_mtDNA],
                                                             ['chr1', 'chrX', 'chrY', 'mtDNA'])]
    # write_vcf(ready_to_map[3])
    threads = min([threads, len(ready_to_map)])
    install_mp_handler()
    pool = mp.Pool(processes=threads)
    pool.map(write_vcf, ready_to_map)
    pool.close()
    pool.join()


def write_pop_file(ts, individuals, populations, target_population, outfile):
    source_populations = [populations[ts.node(ts.individual(indv).nodes[0]).population]
                          if populations[ts.node(ts.individual(indv).nodes[0]).population] != target_population
                          else '-' for indv in individuals]
    with open(outfile, 'w') as f:
        f.write(source_populations[0])
        for pop in source_populations[1:]:
            f.write('\n')
            f.write(pop)
    f.close()


def write_vcf(args):
    """
    Write VCF file for a given chromosome and population
    :param args: (tuple), (TreeSequence, individual IDs, individual names, population, output file, chromosomes)
    """
    ts, individuals, indv_names, outfile, chr, populations = args
    # name the sample
    individual_names = [f"{populations[ts.node(ts.individual(indv).nodes[0]).population]}{indv_name + 1}"
                        for indv, indv_name in zip(individuals, indv_names)]

    # get contig id
    if chr == 'chr1':
        contig_id = '1'
    elif chr == 'chrX':
        contig_id = '23'
    elif chr == 'chrY':
        contig_id = '24'
    elif chr == 'mtDNA':
        contig_id = '26'
    logging.info(f'Writing VCF {outfile}')
    # see https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.write_vcf
    with open(outfile, 'w') as vcf_file:
        ts.write_vcf(vcf_file, individuals=individuals, individual_names=individual_names, contig_id=contig_id)
    logging.info(f"Finished Writing VCF {outfile}")
    vcf_file.close()


def main(argv):
    parser = argparse.ArgumentParser(description="Recaptitate tree sequence and superimpose neutral mutations. "
                                                 "For qualtiy control, recapitation is visualized and SFS is plotted "
                                                 "after superimposing mutations")
    parser.add_argument('--tree_sequences', help='Input tree sequences. Multiple files can be specified separated by '
                                                 'a space. [trees/single_pulse_admixture.trees]',
                        default=['trees/single_pulse_admixture.trees'], nargs='+')
    parser.add_argument('-o', '--output_tree_sequences', help='Output file with mutated tree sequences.'
                                                              'Creates a corresponding pickle file with pointers to '
                                                              'nodes with Y chromosomes and mtDNAs. '
                                                              'Multiple files can be specified separated by a space.'
                                                              '[trees/single_pulse_admixture_mutated.trees]',
                        default=['trees/single_pulse_admixture_mutated.trees'], nargs='+')
    parser.add_argument('--output_vcf_dir', help='Directory where to save VCF file [./data]', default='./data')
    parser.add_argument('-N', '--ancestral_population_size', help='Ancestral population size [7310.370867595234]',
                        type=float, default=7310.370867595234)
    parser.add_argument('--mtDNA_mutation_rate', help='Mutation rate mitochondrial DNA [4e-4]', type=float,
                        default=4e-4)
    parser.add_argument('-m', '--mutation_rate', type=float, default=2.36e-8, help='Germ line mutation rate [2.36e-8')
    parser.add_argument('-l', '--chromosome_length', help='Chromosome Length [1e8]', default=1e8, type=float)
    parser.add_argument('--mtDNA_length', help='Length of mitochondrial DNA [20000]', type=float, default=20000)
    parser.add_argument('--populations', nargs='+', help='Simulated populations in order [afr, eur, ea, amr]',
                        default=['afr', 'eur', 'ea', 'amr'])
    parser.add_argument('-s', '--sample_size', type=int, default=100, help='Sample size, same for each population [100]')
    parser.add_argument('--target_population', help='Population of most interest. Extract full sample size '
                                                    'for this population. [amr]', default='amr')
    parser.add_argument('--off_target_sample_size', help='Reference populations. Use this reduced sample size. [1000]',
                        default=1000, type=int)
    parser.add_argument('--plot_sfs', action='store_true', help='Plot SFS only', default=False)
    parser.add_argument('--recapitate', help='Whether or not do recapitation to ensure complete coalescence [False]',
                        action='store_true', default=False)
    parser.add_argument('--force', help='Force re-run if intermediate files already exist [False]',
                        action='store_true', default=False)
    parser.add_argument('--log', help='Path to log file [log/adding_mutations.log]', default='log/adding_mutations.log')
    parser.add_argument('-t', '--threads', help='Number of processors [16]', default=16, type=int)
    args = parser.parse_args()
    if args.output_vcf_dir.endswith('/'):
        output_vcf_dir = args.output_vcf_dir
    else:
        output_vcf_dir = args.output_vcf_dir + '/'
    if not os.path.isdir(output_vcf_dir):
        os.makedirs(output_vcf_dir)
    logging.basicConfig(filename=args.log, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S')
    # install_mp_handler()
    ancestral_population_size = int(args.ancestral_population_size)
    for ts, ots in zip(args.tree_sequences, args.output_tree_sequences):
        logger = logging.getLogger("Overlaying mutations on {ts}")
        mutated, pointers_Y_mtDNA = add_mutations(ts, ots,
                                                  ancestral_population_size, args.mutation_rate,
                                                  args.mtDNA_mutation_rate, args.chromosome_length, args.mtDNA_length,
                                                  args.populations, args.recapitate, args.force)
        output_base = f"{output_vcf_dir}{ts.split('/')[-1].split('.trees')[0]}"
        extract_population_specific_chromosomes(mutated, pointers_Y_mtDNA, output_base, args.chromosome_length,
                                                args.mtDNA_length, args.populations, args.sample_size,
                                                args.target_population, args.off_target_sample_size, args.plot_sfs,
                                                args.threads)


if __name__ == '__main__':
    main(sys.argv[1:])
