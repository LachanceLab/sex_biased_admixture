#!/usr/bin/env python
import sys
import argparse
import re
import logging
import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import euclidean
from sklearn.cluster import DBSCAN, KMeans, AffinityPropagation, AgglomerativeClustering, SpectralClustering, MeanShift


def loss_clustering(eigenvec, labels, populations):
    """
    Quantify loss based on ancestry of homogeneity of assigned clusters. Distance between individuals of same cluster
    only contributes to loss if they come from different populations
    :param eigenvec: array-like, PCs of samples
    :param labels: array-like, assigned cluster IDs
    :param populations: array-like, populations of origin
    :return: float, loss
    """
    pops = np.array(populations)
    # calculate pairwise distances within clusters
    dist = pairwise_distances(eigenvec)
    # 1 if in same cluster, 0 if in different one
    indicator_matrix_same_cluster =(labels[:, None] == labels[None, :]).astype(int)
    # 0 if sampled from same population, 1 if sampled from different ones
    indicator_matrix_same_pops = (pops[:, None] != pops[None, :]).astype(int)
    # multiple distances by indicator matrices
    loss = np.sum(dist * indicator_matrix_same_pops * indicator_matrix_same_cluster)
    return loss


def find_best_clustering(eigenvec_source, populations_source, not_accepted_algos=[]):
    """
    Find best clustering via GridSearch
    :param eigenvec_source: array-like, PCs of samples of source populations
    :param populations_source: array-like, origin populations of samples
    :param not_accepted_algos: list, clustering algorithms that were previously rejected
    :return: dict, best algorithm + params
    """
    c_error = -1

    # DBSCAN
    if not 'dbscan' in not_accepted_algos:
        logging.info('Running DBSCAN')
        eps_vals = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
        min_samples = np.arange(5, 50, 5)
        for eps in eps_vals:
            for ms in min_samples:
                dbscan = DBSCAN(min_samples=ms, eps=eps).fit(eigenvec_source.values)
                loss = loss_clustering(eigenvec_source, dbscan.labels_, populations_source)
                if loss < c_error or c_error == -1:
                    best_params = {'min_samples': ms, 'eps': eps, 'algo': 'dbscan'}
                    c_error = loss

        logging.info('Best run: Loss: {:.4f}; Params: {}'.format(c_error, best_params))

    n_clusters = np.arange(3, 50, 1)
    if not 'kmeans' in not_accepted_algos:
        # KMeans
        logging.info('Running KMeans')
        for n in n_clusters:
            kmeans = KMeans(n_clusters=n, random_state=1).fit(eigenvec_source.values)
            loss = loss_clustering(eigenvec_source, kmeans.labels_, populations_source)
            if loss < c_error or c_error == -1:
                best_params = {'n_clusters': n, 'algo': 'kmeans'}
                c_error = loss
        logging.info('Best run: Loss: {:.4f}; Params: {}'.format(c_error, best_params))

    if not 'agglomerative' in not_accepted_algos:
        # Agglomerative
        logging.info('Running Agglomerative Clustering')
        for n in n_clusters:
            agglo = AgglomerativeClustering(n_clusters=n).fit(eigenvec_source.values)
            loss = loss_clustering(eigenvec_source, agglo.labels_, populations_source)
            if loss < c_error or c_error == -1:
                best_params = {'n_clusters': n, 'algo': 'agglomerative'}
                c_error = loss
        logging.info('Best run: Loss: {:.4f}; Params: {}'.format(c_error, best_params))

    if not 'affinity' in not_accepted_algos:
        # Affinity
        logging.info('Running Affinity propagation')
        damping = np.arange(0.5, 1.0, 0.05)
        for d in damping:
            affinity = AffinityPropagation(damping=d, max_iter=500, random_state=1).fit(eigenvec_source)
            loss = loss_clustering(eigenvec_source, affinity.labels_, populations_source)
            if loss < c_error  or c_error == -1:
                best_params = {'damping': d, 'algo': 'affinity'}
                c_error = loss

        logging.info('Best run: Loss: {:.4f}; Params: {}'.format(c_error, best_params))

    if not 'meanshift' in not_accepted_algos:
        # MeanShift
        logging.info('Running Mean Shift')
        meanshift = MeanShift().fit(eigenvec_source.values)
        loss = loss_clustering(eigenvec_source, meanshift.labels_, populations_source)
        if loss < c_error or c_error == -1:
            best_params = {'algo': 'meanshift'}
            c_error = loss

        logging.info('Best run: Loss: {:.4f}; Params: {}'.format(c_error, best_params))

    if not 'spectral' in not_accepted_algos:
        # Spectral
        logging.info("Running Spectral Clustering")
        for n in n_clusters:
            spectral = SpectralClustering(n_clusters=n, random_state=1).fit(eigenvec_source.values)
            loss = loss_clustering(eigenvec_source, spectral.labels_, populations_source)
            if loss < c_error or c_error == -1:
                best_params = {'n_clusters': n, 'algo': 'spectral'}
                c_error = loss
        logging.info('Best run: Loss: {:.4f}; Params: {}'.format(c_error, best_params))

    return best_params


def cluster_to_ancestry(labels, populations):
    """
    Map cluster IDs to ancestries based on majority
    :param labels: array-like, cluster IDs
    :param populations: array-like, populations of origin
    :return: dict, mapping cluster id - ancestry
    """
    cluster_to_labels = {}
    for l in labels:
        pops, counts = np.unique(populations[labels == l], return_counts=True)
        c_pop = pops[np.argmax(counts)]
        cluster_to_labels[l] = c_pop
    return cluster_to_labels


def dbscan_predict(dbscan_model, X_test, metric=euclidean):
    """
    Predict DBSCAN cluster belonging
    :param dbscan_model: trained DBSCAN model
    :param X_test: array-like, PCs American samples
    :param metric: sklearn.spatial.distance metric
    :return: array-like, cluster IDs of American samples
    """
    # Taken from https://stackoverflow.com/questions/27822752/scikit-learn-predicting-new-points-with-dbscan
    # Result is noise by default
    y_test = np.ones(shape=len(X_test), dtype=int) * -1

    # Iterate all input samples for a label
    for j, x_new in enumerate(X_test):
        # Find a core sample closer than EPS
        for i, x_core in enumerate(dbscan_model.components_):
            if metric(x_new, x_core) < dbscan_model.eps:
                # Assign label of x_core to x_new
                y_test[j] = dbscan_model.labels_[dbscan_model.core_sample_indices_[i]]
                break
    return y_test


def find_nearest_neighbor(X_train, labels, X_test):
    """
    Define cluster belonging based on nearest neighbor
    :param X_train: array-like, PCs of training samples (non-America)
    :param labels: array-like, cluster ids of training samples
    :param X_test: array-like, PCs of American samples
    :return: array-like, cluster IDs of American samples
    """
    y_test = np.ones(len(X_test))
    for i, x in enumerate(X_test):
        if len(x.shape) == 1:
            x = x[np.newaxis, :]
        # find nearest neighbor
        dist = pairwise_distances(X_train, x)
        y_test[i] = labels[np.argmin(dist)]
    return y_test


def predict_amr_haplogroups(eigenvec_source, eigenvec_amr, populations_source, best_params):
    """
    Predicted cluster belonging of American samples
    :param eigenvec_source: pd.DataFrame, PCs of individuals from source populations
    :param eigenvec_amr: pd.DataFrame, PCs of individuals from American population
    :param populations_source: array-like, population of origin of non-American samples
    :param best_params: dict, best algorithm and parameters
    :return: trained clustering algorithm, predicted cluster ids, corresponding ancestries
    """
    logging.info('Predicting American haplogroups')
    # retrain clustering
    if best_params['algo'] == 'dbscan':
        clustering = DBSCAN(eps=best_params['eps'], min_samples=best_params['min_samples']).fit(eigenvec_source.values)
    elif best_params['algo'] == 'agglomerative':
        clustering = AgglomerativeClustering(n_clusters=best_params['n_clusters']).fit(eigenvec_source.values)
    elif best_params['algo'] == 'kmeans':
        clustering = KMeans(n_clusters=best_params['n_clusters']).fit(eigenvec_source.values)
    elif best_params['algo'] == 'affinity':
        clustering = AffinityPropagation(damping=best_params['damping']).fit(eigenvec_source.values)
    elif best_params['algo'] == 'meanshift':
        clustering = MeanShift().fit(eigenvec_source.values)
    elif best_params['algo'] == 'spectral':
        clustering = SpectralClustering(n_clusters=best_params['n_clusters']).fit(eigenvec_source.values)
    # map cluster labels to ancestry based on majority
    cluster_ancestry_mapping = cluster_to_ancestry(clustering.labels_, populations_source)

    # predict cluster belonging of American individuals
    if best_params['algo'] == 'dbscan':
        amr_cluster = dbscan_predict(clustering, eigenvec_amr.values)
    elif best_params['algo'] == 'agglomerative' or best_params["algo"] == 'spectral':
        amr_cluster = find_nearest_neighbor(eigenvec_source.values, clustering.labels_, eigenvec_amr.values)
    else:
        amr_cluster = clustering.predict(eigenvec_amr.values)
    # convert cluster ids to ancestries
    amr_haplogroups = [cluster_ancestry_mapping[y] for y in amr_cluster]
    # summarize
    ancestry, counts = np.unique([amr_haplogroups], return_counts=True)
    logging.info([f'{counts[i]} {ancestry[i].upper()}' for i in range(len(ancestry))])
    return clustering, amr_cluster, amr_haplogroups


def check_clustering_and_predictions(clustering, predictions, algo):
    """
    Check if clustering and predictions were meaningful
    :param clustering: trained clustering algorithm from sklearn
    :param predictions: np.array, predicted cluster labels for american samples
    :param algo: str, name of algorithm
    :return: boolean, whether or not to accept clustering
    """
    logging.info('Checking clustering and predictions')
    # too many singletons
    if np.unique(clustering.labels_).shape[0] / clustering.labels_.shape[0] > 0.15:
        logging.info(
            '{} found to many clusters. Re-running find_best_clustering with excluding {}'.format(algo, algo))
        return False
    # too many training samples unassigned
    elif clustering.labels_[clustering.labels_ == -1].shape[0] / clustering.labels_.shape[0] > 0.1:
        logging.info('{} left too many samples unassigned. '
                     'Re-running find_best_clustering with excluding {}'.format(algo, algo))
        return False
    # too many american individuals were unassigned
    elif predictions[predictions == -1].shape[0] / predictions.shape[0] > 0.1:
        logging.info('{} left too many American samples unassigned. '
                     'Re-running find_best_clustering with excluding {}'.format(algo, algo))
        return False
    else:
        # accept
        return True


def main(argv):
    parser = argparse.ArgumentParser(description='Cluster ancestral mtDNA and Y haplogroups and assign ancestral '
                                                 'haplogroups of admixed individuals based on proximity to '
                                                 'ancestral clusters')
    parser.add_argument('--eigenvec', help='Path to .eigenvec file generated with plink. Multiple files can be '
                                           'specified separated by a space.', nargs='+')
    parser.add_argument('--log', help='Log file [log/classifying_mtDNA_Y_haplogroups.log]',
                        default='log/classifying_mtDNA_Y_haplogroups.log')
    parser.add_argument('--output_file', help="File to write haplogroup assignments to. Multiple files can be "
                                              "specified separated by a space.", nargs='+')
    args = parser.parse_args()
    logging.basicConfig(filename=args.log, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S')
    logger = logging.getLogger()
    for eigenvec, output_file in zip(args.eigenvec, args.output_file):
        eigenvec = pd.read_csv(eigenvec, sep='\t')
        eigenvec.set_index('#IID', inplace=True)
        populations = np.array([re.sub('[0-9]*', '', ind) for ind in eigenvec.index.values])
        eigenvec_source = eigenvec.loc[populations != 'amr',]
        populations_source = populations[populations != 'amr']
        eigenvec_amr = eigenvec.loc[populations == 'amr',]
        accept_clustering = False
        not_accepted_clustering = []
        while not accept_clustering:
            best_params = find_best_clustering(eigenvec_source, populations_source,
                                               not_accepted_algos=not_accepted_clustering)
            clustering, amr_cluster, amr_haplogroups = predict_amr_haplogroups(eigenvec_source, eigenvec_amr,
                                                                               populations_source, best_params)
            accept_clustering = check_clustering_and_predictions(clustering, amr_cluster, best_params['algo'])
            if not accept_clustering:
                logging.info('Did not accept clustering {}'.format(best_params))
                not_accepted_clustering.append(best_params['algo'])
            if len(not_accepted_clustering) == 6:
                logging.info('None of the tried clustering algorithms was accepted ({}). '
                             'Try something else.'.format(not_accepted_clustering))
                sys.exit(1)
        logging.info('Accepted {}'.format(best_params))
        cluster_ancestry_mapping = cluster_to_ancestry(clustering.labels_, populations_source)
        source_haplogroups = [cluster_ancestry_mapping[y] for y in clustering.labels_]
        haplogroups = pd.DataFrame(index=eigenvec.index.values)
        haplogroups.loc[eigenvec_source.index.values, 'haplogroup'] = source_haplogroups
        haplogroups.loc[eigenvec_amr.index.values, 'haplogroup'] = amr_haplogroups
        logging.info(f'Writing haplogroup assignments to {output_file}')
        haplogroups.to_csv(output_file, sep='\t', header=True, index=True)


if __name__ == '__main__':
    main(sys.argv[1:])
