#! /usr/bin/env python3

"""Ad hoc script: get coverage data
"""

import sys
import os
import glob
import re
import shlex
import argparse
import logging
import subprocess
import csv
import json
import pandas as pd


logger = logging.getLogger(__name__)

        
        
def get_args(args):
    
    from collections import defaultdict
    
    final_options = defaultdict(list)
    final_options_str = ''
    if args is not None:
        for opt in args:
            final_options['{}'.format(opt.rsplit()[0].strip('--'))].append(' '.join(opt.rsplit()[1:]))

        final_options_str = ' '.join(['--{} {}'.format(key, val) for (key, values) in
                                      final_options.items() for val in values])
    return final_options_str


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Run bamCoverage on a set of files.""")

    parser.add_argument('flist', help="The list of files, with full path in the first column.")

    parser.add_argument('output', help="""The output directory.""")
    
    parser.add_argument('--dt-options', help="""A space-delimited
        list of options to pass to bamCoverage. Each option
        must be quoted separately as in "--option value", using soft 
        quotes, where '--option' is the long parameter name, and 'value'
        is the value given to this parameter.""", nargs='*', type=str)
    
    parser.add_argument('--overwrite', help='''If this flag is present, then existing files
        will be overwritten.''', action='store_true')

    add_logging_options(parser)
    args = parser.parse_args()
    update_logging(args)

    # check output path
    if os.path.exists(args.output):
        args.output = os.path.join(args.output, '')
    else:
        msg = "Invalid output path or wrong permission: {}. Terminating.".format(args.output)
        raise OSError(msg)
    
    msg = "[get_coverage]: {}".format(' '.join(sys.argv))
    logger.info(msg)
    
    # get bamCoverage arguments, no check on validity of arguments
    bamc_args = get_args(args.dt_options)
    
    file_list = pd.read_csv(args.flist,
                            header=None,
                            usecols=[0],
                            names=['path'])
    
    for _, row in file_list.iterrows():
        
        output_filename = '{}.bw'.format(os.path.splitext(os.path.basename(row.path))[0])
        output_filename = os.path.join(args.output, output_filename)
        cmd = 'bamCoverage -b {} -o {} -of bigwig --ignoreDuplicates {}'.format(row.path,
                                                                                output_filename,
                                                                                bamc_args) 
        #cmd = shlex.split(cmd)
        # check if the output files exist
        out_files = [output_filename]
        all_out_exists = all([os.path.exists(of) for of in out_files])
        if args.overwrite or not all_out_exists:
            msg = 'Processing {}'.format(row.path)
            logger.info(msg)
            #p = subprocess.Popen(cmd) 
            ret_code = subprocess.call(cmd, shell=True)
            if ret_code != 0:
                raise subprocess.CalledProcessError(ret_code, cmd)


if __name__ == '__main__':
    main()
