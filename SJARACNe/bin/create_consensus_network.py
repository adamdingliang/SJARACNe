#!/usr/bin/env python3

import sys
import os
import argparse
import math
import logging
import numbers
import numpy as np
import pathlib
import re
from scipy import stats
import pandas as pd


DEFAULT_MIN_RECURRENCE = 6
_UNSET = object()


def minimum_recurrence(value):
    """Parse a positive integer recurrence threshold for argparse."""
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError('minimum recurrence must be an integer')
    if parsed < 1:
        raise argparse.ArgumentTypeError('minimum recurrence must be at least 1')
    return parsed


def consensus_probability(value):
    """Parse the deprecated legacy consensus probability threshold."""
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError('consensus p-value must be a number')
    if not math.isfinite(parsed) or not 0.0 < parsed <= 1.0:
        raise argparse.ArgumentTypeError('consensus p-value must be within (0, 1]')
    return parsed


def main():
    """ Handles arguments and invokes the driver function. """
    head_description = '''Create a consensus network based on SJARACNe bootstrap networks.'''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=head_description)
    parser.add_argument('-a', '--adjmat-dir', metavar='STR', required=True, help='directory with adjacent matrix')
    selection_group = parser.add_mutually_exclusive_group()
    selection_group.add_argument(
        '-k', '--min-recurrence', metavar='INT', type=minimum_recurrence,
        help=(
            'minimum number of distinct bootstrap networks containing an edge; '
            'default: 6'
        ),
    )
    selection_group.add_argument(
        '-p', '--p-value', metavar='FLOAT', type=consensus_probability,
        help=(
            'deprecated legacy normal-approximation p-value threshold; use '
            '--min-recurrence instead'
        ),
    )
    parser.add_argument('-e', '--exp-mat', metavar='STR', required=True, help='expression matrix file')
    parser.add_argument('-o', '--out-dir', metavar='STR', required=True, help='output directory')
    parser.add_argument('-s', '--subnet', metavar='STR', help='file with gene symbols of interest to build a subnet')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    logging.basicConfig(level=logging.INFO)
    logging.info('Create an initial consensus network ...')
    selection = {}
    if args.min_recurrence is not None:
        selection['min_recurrence'] = args.min_recurrence
    elif args.p_value is not None:
        logging.warning(
            'Legacy consensus p-value filtering uses a normal approximation; '
            'prefer --min-recurrence.'
        )
        selection['p_value'] = args.p_value
    network = create_consensus_network(
        args.adjmat_dir,
        out_dir=args.out_dir,
        **selection,
    )
    logging.info('Done')
    logging.info('Create an enhanced consensus network ...')
    create_enhanced_consensus_network(args.exp_mat, network, args.out_dir, args.subnet)
    logging.info('All done')


def create_consensus_network(
    adjmat_dir,
    p_value=_UNSET,
    out_dir=None,
    *,
    min_recurrence=None,
):
    """ Create a consensus network based on SJARACNe bootstrap networks
    Args:
        adjmat_dir: directory with adjacent matrix
        p_value: deprecated legacy normal-approximation p-value threshold. An
            explicit None retains the historical Bonferroni behavior.
        out_dir: output directory
        min_recurrence: minimum number of distinct bootstrap networks that
            must contain an ordered edge. When neither selection mode is
            supplied, defaults to 6.
    Returns:
        none
    """
    legacy_probability_mode = p_value is not _UNSET
    if legacy_probability_mode and min_recurrence is not None:
        raise ValueError('p_value and min_recurrence are mutually exclusive')

    if not legacy_probability_mode and min_recurrence is None:
        min_recurrence = DEFAULT_MIN_RECURRENCE

    if min_recurrence is not None:
        if (
            isinstance(min_recurrence, bool)
            or not isinstance(min_recurrence, numbers.Integral)
        ):
            raise ValueError('minimum recurrence must be an integer')
        min_recurrence = int(min_recurrence)
        if min_recurrence < 1:
            raise ValueError('minimum recurrence must be at least 1')

    if legacy_probability_mode and p_value is not None:
        try:
            parsed_p_value = float(p_value)
        except (TypeError, ValueError):
            raise ValueError('consensus p-value must be a number')
        if not math.isfinite(parsed_p_value) or not 0.0 < parsed_p_value <= 1.0:
            raise ValueError('consensus p-value must be within (0, 1]')
        p_value = parsed_p_value

    if out_dir is None:
        raise ValueError('out_dir is required')

    total_edge_in_runs = []
    bootstrap_run_num = 0
    parameters = ''
    total_edge_number = {}
    total_mi = {}

    # Processing all bootstrap networks, summarizing them into corresponding variables
    for adj_file in sorted(os.listdir(adjmat_dir)):
        edges_in_run = {}
        duplicate_edge_count = 0
        # Opening each bootstrap file
        with open(pathlib.PurePath(adjmat_dir).joinpath(adj_file), "r") as fadj:
            for line in fadj:
                # Processing header lines
                if line[0] == '>' and bootstrap_run_num == 0:
                    parameters += line
                # Processing non header lines representing the network
                if line[0] != '>':
                    tokens = line.split('\t')
                    # Tokenizing each non header line with tab delimiter
                    hub_id = tokens[0]  # First token is the hub id
                    # Iterating on all adjacent genes: Odd indexes are the connected genes and even indexes are
                    # the corresponding value to the edge between the hub gene and the gene with an odd index
                    # appearing before the value in the tokens list
                    for index in range(1, len(tokens), 2):
                        key = (hub_id, tokens[index])  # Ordered source-target edge
                        mi = float(tokens[index + 1])

                        if not math.isfinite(mi):
                            raise ValueError(
                                "Non-finite MI value for edge {}----{} in bootstrap "
                                "file '{}': {}".format(key[0], key[1], adj_file, mi)
                            )

                        if key in edges_in_run:
                            duplicate_edge_count += 1
                            if edges_in_run[key] != mi:
                                raise ValueError(
                                    "Conflicting MI values for duplicate edge {}----{} in bootstrap "
                                    "file '{}': {} versus {}".format(
                                        key[0], key[1], adj_file, edges_in_run[key], mi
                                    )
                                )
                            continue

                        edges_in_run[key] = mi

        if duplicate_edge_count:
            logging.warning(
                "Ignored %d duplicate edge occurrence(s) in bootstrap file '%s'",
                duplicate_edge_count,
                adj_file,
            )

        total_edge_in_runs.append(len(edges_in_run))
        for key, mi in edges_in_run.items():
            # A bootstrap is one Bernoulli observation for an edge, regardless
            # of how many times a malformed file repeats that edge.
            total_edge_number[key] = total_edge_number.get(key, 0) + 1
            total_mi[key] = total_mi.get(key, 0.0) + mi

        # Increment the bootstrap file index
        bootstrap_run_num += 1

    if bootstrap_run_num == 0:
        raise ValueError(
            "No bootstrap adjacency files found in '{}'.".format(adjmat_dir)
        )

    if min_recurrence is not None and min_recurrence > bootstrap_run_num:
        raise ValueError(
            'minimum recurrence {} exceeds the number of bootstrap networks {}'
            .format(min_recurrence, bootstrap_run_num)
        )

    # Do not leave an output directory behind when input or recurrence
    # validation fails.
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    edge_count = len(total_edge_number)
    p_threshold = None
    mu = None
    sigma = None

    if legacy_probability_mode:
        mu = 0.0
        variance = 0.0
        # Preserve the historical normal approximation for explicit legacy use.
        if edge_count > 0:
            for edge_total in total_edge_in_runs:
                prob = float(edge_total) / float(edge_count)
                mu += prob
                variance += prob * (1 - prob)
        sigma = np.sqrt(variance)

        if edge_count > 0:
            bonferroni_alpha = 0.05 / edge_count
            bonferroni_alpha_text = str(bonferroni_alpha)
        else:
            bonferroni_alpha = None
            bonferroni_alpha_text = 'N/A (no edges tested)'

        p_threshold = bonferroni_alpha if p_value is None else p_value

        with open(pathlib.PurePath(out_dir).joinpath('bootstrap_info_.txt'), 'w') as f_info:
            f_info.write('Total edge tested: {}\n'.format(str(edge_count)))
            f_info.write(
                'Bonferroni corrected (0.05) alpha: {}\n'.format(
                    bonferroni_alpha_text
                )
            )
            f_info.write('mu: {}\n'.format(str(mu)))
            f_info.write('sigma: {}\n'.format(str(sigma)))
    else:
        recurrence_fraction = float(min_recurrence) / float(bootstrap_run_num)
        with open(pathlib.PurePath(out_dir).joinpath('bootstrap_info_.txt'), 'w') as f_info:
            f_info.write('Total edge tested: {}\n'.format(str(edge_count)))
            f_info.write('Consensus selection: minimum recurrence\n')
            f_info.write('Bootstrap networks: {}\n'.format(bootstrap_run_num))
            f_info.write('Minimum recurrence: {}\n'.format(min_recurrence))
            f_info.write(
                'Minimum recurrence fraction: {}\n'.format(
                    '{:.12g}'.format(recurrence_fraction)
                )
            )

    # Writing out the parameters that the bootstrap networks are constructed with plus other
    # parameters that is used to create consensus network
    parameters += '>  Bootstrap No: {}\n'.format(str(bootstrap_run_num))
    if min_recurrence is not None:
        parameters += '>  Consensus selection minimum recurrence\n'
        parameters += '>  Minimum recurrence {} of {}\n'.format(
            min_recurrence,
            bootstrap_run_num,
        )
        parameters += '>  Minimum recurrence fraction {}\n'.format(
            '{:.12g}'.format(float(min_recurrence) / float(bootstrap_run_num))
        )
    parameters += '>  Source: sjaracne2\n'
    out_network_path = pathlib.PurePath(out_dir).joinpath('consensus_network_3col_.txt')
    parameters += '>  Output network: {}\n'.format(out_network_path)
    with open(pathlib.PurePath(out_dir).joinpath('parameter_info_.txt'), 'w') as parameter_file:
        parameter_file.write(parameters)

    # Write the consensus network using the selected recurrence or legacy
    # probability rule.
    with open(out_network_path, 'w') as f_consensus_network:
        header = 'source\ttarget\tMI\n'
        f_consensus_network.write(header)

        # Iterate over all edges in a sorted fashion
        for key in sorted(total_edge_number.keys()):
            # Extract first two gene involving an edge from the key (edge)
            gene1, gene2 = key

            if min_recurrence is not None:
                retain_edge = total_edge_number[key] >= min_recurrence
            else:
                # Explicit legacy mode: retain the historical normal-tail gate.
                z = (
                    float(total_edge_number[key] - mu) / float(sigma)
                    if sigma != 0
                    else 100
                )
                retain_edge = uprob(z) < p_threshold

            if retain_edge:
                # Computing MI corresponding to an edge remaining in the network
                mi = '{0:.4f}'.format(float(total_mi[key]) / float(total_edge_number[key]))
                f_consensus_network.write('{}\t{}\t{}\n'.format(gene1, gene2, mi))
    return out_network_path


def create_enhanced_consensus_network(exp_mat, network, out_dir, subnet=None):
    """ Add more information to a consensus network to create an enhanced network.
    Args:
        exp_mat (str): path to an expression matrix file
        network (str): path to a consensus network file
        out_dir (str): path to an output directory
        subnet (str, optional): path to a gene symbol file
    Return:
        None
    """
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # Build output file name based on input network file path
    network_file_name = os.path.basename(network)
    input_net_name_tokens = re.split(r"_|\.", network_file_name)
    out_file_name = ("_".join(input_net_name_tokens[0:2]) + "_ncol_" +
                     "_".join(input_net_name_tokens[3:len(input_net_name_tokens) - 1]))
    header = ("source", "target", "source.symbol", "target.symbol", "MI", "pearson",
              "spearman", "slope", "p-value")

    # Read in subnetwork file
    subnet_list = []
    out_subnet = None
    if subnet:
        out_subnet = open(pathlib.PurePath(out_dir).joinpath("sub_" + out_file_name + ".txt"), "w")
        with open(subnet) as subnet_file:
            for _id in subnet_file:
                subnet_list.append(_id.split("\n")[0].strip())

    exp = pd.read_csv(exp_mat, sep="\t", index_col=0)

    exp_dict = dict()
    with open(network, 'r') as fnet:
        fnet.readline()
        with open(pathlib.PurePath(out_dir).joinpath(out_file_name + ".txt"), 'w') as fout:
            fout.write('\t'.join(header) + '\n')
            for line in fnet:
                tokens = line.split('\t')
                node1 = tokens[0]
                node2 = tokens[1]
                mi = float(tokens[2])

                if node1 in exp_dict:
                    temp = exp_dict[node1]
                    gene_symbol1 = temp[0]
                    exp_values1 = temp[1]
                else:
                    exp_symbol_values1 = exp.loc[node1].values
                    gene_symbol1 = str(exp_symbol_values1[0])
                    exp_values1 = exp_symbol_values1[1:].astype(float)
                    exp_dict[node1] = (gene_symbol1, exp_values1)

                if node2 in exp_dict:
                    temp = exp_dict[node2]
                    gene_symbol2 = temp[0]
                    exp_values2 = temp[1]
                else:
                    exp_symbol_values2 = exp.loc[node2].values
                    gene_symbol2 = str(exp_symbol_values2[0])
                    exp_values2 = exp_symbol_values2[1:].astype(float)
                    exp_dict[node2] = (gene_symbol2, exp_values2)

                slope, intercept, r, p, stderr = stats.linregress(exp_values1, exp_values2)
                scc, sp = stats.spearmanr(exp_values1, exp_values2)
                pcc, pp = stats.pearsonr(exp_values1, exp_values2)

                row = (node1, node2, gene_symbol1, gene_symbol2,
                       "{0:.4f}".format(mi), "{0:.4f}".format(pcc), "{0:.4f}".format(scc),
                       "{0:.4f}".format(slope), "{0:.4f}".format(p))
                fout.write('\t'.join(row) + '\n')

                if out_subnet is not None:
                    if gene_symbol1 in subnet_list or gene_symbol2 in subnet_list:
                        out_subnet.write('\t'.join(row) + '\n')

    if out_subnet is not None:
        out_subnet.close()


def uprob(n):
    """ Implemented in statistics.py module inspired by Statistics::Distributions::uprob function in perl.
    Args:
        n (float): z-score
    Returns:
        p (float): p value
    """
    p = 0
    if abs(n) < 1.9:
        p = (1 + abs(n) * (0.049867347 + abs(n) * (
            0.0211410061 + abs(n) * 0.0032776263 + abs(n)
            * (0.0000380036 + abs(n) * (0.0000488906 + abs(n) * 0.000005383))))) ** (-16) / 2
    elif abs(n) <= 100:
        for i in range(18, 0, -1):
            p = i / (abs(n) + p)
        p = math.exp(-0.5 * abs(n) * abs(n)) / math.sqrt(2 * math.pi) / (abs(n) + p)
    if n < 0:
        p = 1 - p
    return p


if __name__ == '__main__':
    main()
