#! /usr/bin/env python3

"""Convert BAM files to bigWig by running deepTools bamCoverage.

Functions:
    _get_args
"""

import sys
import os
import glob
import argparse
import logging
import pandas as pd

import phub.utils as utils

logger = logging.getLogger(__name__)

        
def _get_args(args):
    
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

    parser.add_argument('inputDir', help="""The input directory with BAM files, or path to
        a file (txt or yaml) containing full path for each BAM file. The input format is 
        specified with [--input-format/-fmt]. All BAM files MUST have the extension bam
        [--input-format glob], and names are extracted based on this information, unless
        a yaml file is used, in which case the key:value pair is used to assign names.""")
    
    parser.add_argument('outputDir', help="""The output directory. BigWig files are 
        written to this location.""")
    
    
    parser.add_argument('-fmt', '--input-format', help="""The 'type' of input, either 'glob', 
        'txt', or 'yaml' file. With 'glob', regex pattern can be specified with [--pattern].
        With 'yaml', the key must be specified with [--key].""", type=str, 
        choices=['glob', 'txt', 'yaml'], default='glob')

    parser.add_argument('--pattern', help="""A comma-separated list of patterns used to glob
        BAM files.""", default='', type=str, nargs='*')
    
    parser.add_argument('--key', help="""The key to access the list of files if using
        a yaml input file.""", default='samples', type=str)
    
    # TODO: some arguments won't be possible to pass, include defaults?
    parser.add_argument('--dt-options', help="""A space-delimited list of options to pass 
        to bamCoverage (deepTools). Each option must be quoted separately as in 
        "--option value", using soft quotes, where '--option' is the long parameter name, 
        and 'value' is the value given to this parameter. There is no check on the validity
        of arguments.""", nargs='*', type=str)
  
    utils.add_file_options(parser)
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)
  
    msg = "[get_bam2bigWig]: {}".format(' '.join(sys.argv))
    logger.info(msg)
    
    # check that deepTools bamCoverage is callable
    programs = ['bamCoverage']
    utils.check_programs_exist(programs)

    # check output path
    if os.path.exists(args.outputDir):
        args.outputDir = os.path.join(args.outputDir, '')
    else:
        msg = "Invalid output path or wrong permission: {}. Terminating.".format(args.outputDir)
        raise OSError(msg)
    
    # get bamCoverage arguments
    bamc_args = _get_args(args.dt_options)
    
    # input format - fetch BAM files
    if args.input_format == 'glob': 
        filenames = list(_glob_re(r'.*(bam)', os.listdir(args.inputDir)))
        match = r'.*({}).*'.format('|'.join(args.pattern))
        filenames = list(_glob_re(match, filenames))
        filenames = [os.path.join(args.inputDir, f) for f in filenames]
    elif args.input_format == 'txt':
        filenames = pd.read_csv(args.inputDir,
                                header=None,
                                usecols=[0],
                                names=['f']).f.to_list()
    elif args.input_format == 'yaml':
        config = yaml.load(open(args.inputDir), Loader=yaml.FullLoader)
        try:
            filenames = list(config[args.key].values())
            fileMapping = config[args.key]
        except:
            msg = 'Missing/wrong key: [--input-format yaml] but [--key] is wrong!'
            raise KeyError(msg)
            
    if not args.input_format == 'yaml':
        fileMapping = {}
        for f in filenames:
            # TODO: interaction POSIX - Windows styled path...?
            outputFile = os.path.splitext(os.path.basename(f))[0]
            fileMapping[outputFile] = f
    
    for bigWig, bam in fileMapping.items():
        if not os.path.exists(bam):
            msg = "Could not find the BAM file: {}. Terminating.".format(bam)
            raise FileNotFoundError(msg)        
        
        filename = os.path.join(args.outputDir, '{}.bw'.format(bigWig))
        in_files = [bam]
        out_files = [filename]
        cmd = 'bamCoverage -b {} -o {} -of bigwig {}'.format(bam,
                                                             filename,
                                                             bamc_args) 
        utils.call_if_not_exists(cmd,
                                 out_files,
                                 in_files=in_files,
                                 overwrite=args.overwrite,
                                 call=True)
    

if __name__ == '__main__':
    main()
