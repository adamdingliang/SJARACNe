#!/usr/bin/env python3

import os
import sys
import argparse
import json
import subprocess
import shlex
import logging
import pathlib


DEFAULT_MIN_RECURRENCE = 6


def sampling_fraction(value):
    parsed = float(value)
    if not 0.0 < parsed <= 1.0:
        raise argparse.ArgumentTypeError('subsample fraction must be within (0, 1]')
    return parsed


def sampling_size(value):
    parsed = int(value)
    if parsed < 2:
        raise argparse.ArgumentTypeError('subsample size must be at least 2')
    return parsed


def minimum_recurrence(value):
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError('minimum recurrence must be an integer')
    if parsed < 1:
        raise argparse.ArgumentTypeError('minimum recurrence must be at least 1')
    return parsed


def bootstrap_count(value):
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError('bootstrap count must be an integer')
    if parsed < 1:
        raise argparse.ArgumentTypeError('bootstrap count must be at least 1')
    return parsed


def consensus_probability(value):
    parsed = float(value)
    if not 0.0 < parsed <= 1.0:
        raise argparse.ArgumentTypeError('consensus p-value must be within (0, 1]')
    return parsed


def main():
    head_description = '''SJARACNe is a scalable tool for gene network reverse engineering.'''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=head_description)

    # Create a parent parser with common arguments for every subparser
    parent_parser = argparse.ArgumentParser(description='hello', add_help=False)
    parent_parser.add_argument('-e', '--exp-file', metavar='FILE', required=True,
                               help='Path to an expression matrix file, row indexes are used as the'
                               'nodes in the network.')
    parent_parser.add_argument('-g', '--hub-genes', metavar='FILE', required=True,
                               help='Path to a file containing a list of symbols to be considered as hub genes.')
    consensus_group = parent_parser.add_mutually_exclusive_group()
    consensus_group.add_argument(
        '-k', '--min-recurrence', metavar='INT', type=minimum_recurrence,
        help=(
            'Minimum number of distinct resampled networks containing an edge; '
            'default: 6.'
        ),
    )
    consensus_group.add_argument(
        '-pc', '--p-value-consensus', metavar='FLOAT', type=consensus_probability,
        help=(
            'Deprecated legacy normal-approximation P-value threshold for the '
            'consensus network.'
        ),
    )
    parent_parser.add_argument('-pb', '--p-value-bootstrap', metavar='FLOAT', default=1e-7,
                               help='P-value threshold to filter mutual information in each resampled network.')
    parent_parser.add_argument(
        '-M', '--apmi-null-model', metavar='FILE',
        help='Exact-m, exact-depth estimator-matched AP-MI GPD-tail model passed to every resampled network.',
    )
    parent_parser.add_argument('-d', '--depth', metavar='INT', default=40, help='maximum partitioning depth.')
    parent_parser.add_argument('-c', '--config-dir', metavar='DIR', help='Directory containing ARACNe configuration '
                                                                         'files. Use default configs if not provided.')
    parent_parser.add_argument('-n', '--bootstrap-num', metavar='INT', type=bootstrap_count, default=100,
                               help='Number of resampled networks to generate (legacy option name).')
    sampling_group = parent_parser.add_mutually_exclusive_group()
    sampling_group.add_argument(
        '-sf', '--subsample-fraction', metavar='FLOAT', type=sampling_fraction,
        help='Fraction of eligible observations sampled without replacement in each network; default: 0.8.',
    )
    sampling_group.add_argument(
        '-sm', '--subsample-size', metavar='INT', type=sampling_size,
        help='Exact number of eligible observations sampled without replacement in each network.',
    )
    parent_parser.add_argument('-o', '--output-dir', metavar='DIR', required=True,
                               help='Path to final output directory.')
    parent_parser.add_argument('-tmp', '--tmpdir-prefix', dest='tmpdir_prefix',metavar='DIR', required=True,
                               help='Specify tmp path,default is /tmp.')

    subparsers = parser.add_subparsers(title='Subcommands', help='platforms', dest='subcommand')
    subparsers.required = True
    # Create a subparser for running cwltool
    subparser_local = subparsers.add_parser('local', parents=[parent_parser], help='run cwltool in a local workstation')
    subparser_local.add_argument('-s', '--serial', help='run cwltool in serial mode', action='store_true')

    # Create a subparser for running cwlexec
    subparser_lsf = subparsers.add_parser('lsf', parents=[parent_parser], help='run cwlexec in a IBM LSF cluster')
    subparser_lsf.add_argument('-j', '--config-json', metavar='FILE', required=True, help='LSF-specific configuration '
                                                                                          'file in JSON format to be '
                                                                                          'used for workflow execution')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    recurrence = None
    if args.p_value_consensus is None:
        recurrence = (
            DEFAULT_MIN_RECURRENCE
            if args.min_recurrence is None
            else args.min_recurrence
        )
        if recurrence > args.bootstrap_num:
            parser.error(
                'minimum recurrence {} exceeds bootstrap count {}'.format(
                    recurrence,
                    args.bootstrap_num,
                )
            )

    # to make executable and config findable
    installed_path = os.path.dirname(os.path.realpath(__file__))
    os.environ['PATH'] += (os.pathsep + installed_path + '/bin')
    cwl_path = installed_path + '/cwl'
    default_config_path = installed_path + '/config'
    if args.config_dir:
        config_dir = args.config_dir
    else:
        config_dir = default_config_path

    if args.subsample_size is not None:
        subsample_spec = str(args.subsample_size)
    else:
        fraction = 0.8 if args.subsample_fraction is None else args.subsample_fraction
        subsample_spec = '{:.12g}%'.format(fraction * 100.0)
        
    if not os.path.isdir(args.tmpdir_prefix):
        os.makedirs(args.tmpdir_prefix)    

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    output_dir_name = os.path.basename(args.output_dir)
    if args.p_value_consensus is not None:
        consensus_selection = 'p_value_consensus: {}\n'.format(
            args.p_value_consensus
        )
    else:
        consensus_selection = 'min_recurrence: {}\n'.format(recurrence)
    # Create input yml file in a temp directory
    with open(pathlib.PurePath(args.output_dir).joinpath('sjaracne_workflow.yml'), 'w') as fp_yml:
        logging.info(fp_yml.name)
        contents = 'exp_file:\n  class: File\n  path: {}\n' \
                   'probe_file:\n  class: File\n  path: {}\n' \
                   '{}' \
                   'p_value_bootstrap: {}\n' \
                   'depth: {}\n' \
                   'aracne_config_dir:\n  class: Directory\n  path: {}\n' \
                   'bootstrap_num: {}\n' \
                   'subsample_spec: {}\n' \
                   'final_out_dir_name: {}'.format(os.path.abspath(args.exp_file), os.path.abspath(args.hub_genes),
                                                   consensus_selection, args.p_value_bootstrap, args.depth,
                                                   config_dir, args.bootstrap_num, json.dumps(subsample_spec),
                                                   output_dir_name)
        if args.apmi_null_model is not None:
            contents += '\napmi_null_model:\n  class: File\n  path: {}'.format(
                json.dumps(os.path.abspath(args.apmi_null_model))
            )
        logging.info(contents)
        fp_yml.write(contents)
        fp_yml.flush()
        fp_yml.seek(0)

        if args.subcommand == 'local':
            if args.serial:
                cmd = 'cwltool --tmpdir-prefix {} --outdir {} {}/sjaracne_workflow.cwl {}'.format(args.tmpdir_prefix, args.output_dir, cwl_path,fp_yml.name)
            else:
                cmd = 'cwltool --tmpdir-prefix {} --parallel --outdir {} {}/sjaracne_workflow.cwl {}'.format(args.tmpdir_prefix, args.output_dir,cwl_path, fp_yml.name)
        elif args.subcommand == 'lsf':
                cmd = 'cwlexec -pe PATH -c {} --outdir {} {}/sjaracne_workflow.cwl {}'.format(
                    args.config_json, args.output_dir, cwl_path, fp_yml.name)
        else:
            sys.exit('Error - invalid subcommand.')
        logging.info(cmd)
        run_shell_command_call(cmd)

    logging.info('All done.')


def run_shell_command_call(cmd):
    """ Wrapper of subprocess.check_call to take a cmd string as input
    Args:
        cmd (str): command to run
    """
    cmd_to_exec = shlex.split(cmd)
    subprocess.check_call(cmd_to_exec)


if __name__ == "__main__":
    main()
