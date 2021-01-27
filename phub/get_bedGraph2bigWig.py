#! /usr/bin/env python3

"""Convert BED4+/bedGraph files to bigWig.

Note: Requires the executable "bedGraphToBigWig" and the "chrom.sizes"
      file which can be obtained with the "fetchChromSizes" script.
      The chrom names from the predictions must match those
      of chrom.sizes (UCSC), otherwise they must be re-written by
      passing the right options.

Functions:
    _glob_re
    _read_bed
    _get_bed
    _convert
"""

import sys
import os
import glob
import re
import argparse
import logging
import yaml
import csv
import json

import pandas as pd

import phub.utils as utils

logger = logging.getLogger(__name__)


default_fields = [
    "chrom", "chromStart", "chromEnd", "name", "score", "strand",
    "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "chromStarts"
]
BEDGRAPH_FIELDS = default_fields[0:3]
BEDGRAPH_FIELDS.extend(['dataValue'])


def _glob_re(pattern, strings):
    return filter(re.compile(pattern).match, strings)


def _read_bed(filename, header, args, sep='\t', **kwargs):
    
    """Reads a BED/bedGraph file into a pandas data frame. This function assumes that 
    field names are prepended with a comment character (see pbio). Otherwise
    a header is added as specified.
    """
    
    if not args.no_header:
        bed = pd.read_csv(filename, sep=sep, **kwargs)
        bed.columns = [c.replace("#", "") for c in bed.columns]
    else:
        bed = pd.read_csv(filename, sep=sep, header=None, **kwargs)
        num_columns = len(bed.columns)
        bed.columns = header[:num_columns]

    # either way, make sure the first column (chrom) is treated as a string
    chrom_name = bed.columns[0]
    chrom_column = bed[chrom_name]
    bed[chrom_name] = chrom_column.astype(str)

    return bed


def _get_bed(filename, output_filename, args, strand=None):

    """Get BED/bedGraph file and adjust features.
    """

    bed = _read_bed(filename, BEDGRAPH_FIELDS, args)
    
    # get fields
    header = list(bed.columns)
    chromField = header[0]
    startField = header[1]
    endField = header[2]
    dataValueField = header[3]
    strandField = header[5]
    if args.value_field is not None:
        dataValueField = args.value_field
    fields_to_keep = [chromField, startField, endField, dataValueField]
    
    # split if strand is given
    if strand is not None:
        bed = bed[bed[strandField]==strand]
              
    # adjust chrom field
    if args.add_chr:
        bed[chromField] = 'chr' + bed[chromField].astype(str)
    if args.chr_dict:
        for chrom_old, chrom_new in args.chr_dict.items():
            seqname_m = bed[chromField] == str(chrom_old)
            bed.loc[seqname_m, chromField] = str(chrom_new)
    if args.chr_file:
        chr_map = pd.read_csv(args.chr_file, 
                              header=None,
                              index_col=0, 
                              squeeze=True).to_dict()
        bed.replace({chromField: chr_map}, inplace=True)

    # sort on the chrom field, and then on the chromStart field.
    bed.sort_values([chromField, startField], ascending=[True, True], inplace=True)

    # remove unused fields
    bed = bed[fields_to_keep]

    # write bedGraph file to output directory
    tmp = os.path.join(args.outputDir, '{}.bedGraph'.format(output_filename))
    if os.path.exists(tmp) and not args.overwrite:
        msg = "Temporary output file {} already exists. Skipping.".format(tmp)
        logger.warning(msg)
    else:
        bed.to_csv(tmp,
                   sep='\t',
                   index=False,
                   header=False,
                   quoting=csv.QUOTE_NONE)
    
    return tmp

        
def _convert(bg, bw, args):

    in_files = [bg, args.chrSizes]
    out_files = [bw]
    cmd = "bedGraphToBigWig {} {} {}".format(bg,
                                             args.chrSizes,
                                             bw)
    utils.call_if_not_exists(cmd,
                             out_files,
                             in_files=in_files,
                             overwrite=args.overwrite,
                             call=True)
    if not args.keep:
        try:
            os.remove(bg)
            msg = "Removing: {}".format(bg)
            logger.info(msg)
        except OSError:
            msg = "Could not remove: {}".format(bg)
            logger.info(msg)
            
            
def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Convert bedGraph to bigWig files by 
        calling the executable program 'bedGraphToBigWig'. The executable must be available 
        on the user's path.""")

    parser.add_argument('inputDir', help="""The input directory with bedGraph files, or path 
        to a file (txt or yaml) containing full path for each file to convert. The input
        format is specified with [--input-format/-fmt]. The bedGraph extension is specified
        by [--extension], and names are extracted based on this information, unless a 
        yaml file is used, in which case the key:value pair is used to assign names.""")
    
    parser.add_argument('outputDir', help="""The output directory. All bedGraph files are 
        also temporarily re-written to this location.""")

    parser.add_argument('chrSizes', help="The 'chrom.sizes' file for the UCSC database.")

    
    parser.add_argument('-fmt', '--input-format', help="""The 'type' of input, either 'glob', 
        'txt', or 'yaml' file. With 'glob', regex pattern can be specified with [--pattern].
        With 'yaml', the key must be specified with [--key].""", type=str, 
        choices=['glob', 'txt', 'yaml'], default='glob')

    parser.add_argument('--pattern', help="""A comma-separated list of patterns used to glob
        BED files.""", default='', type=str, nargs='*')
    
    parser.add_argument('--key', help="""The key to access the list of files if using
        a yaml input file.""", default='samples', type=str)
    
    parser.add_argument('--extension', help="""The input bedGraph extension before any 
        compression (gz).""", type=str, default='bedGraph')
    
    parser.add_argument('--skip', help="""If this flag is present then temporary bedGraph
        files are not generated, and the input files are converted directly. There is no format
        check! Options/flags such as [--add-chr, --chr-dict, --chr-file, --no-header, 
        --split-strand, --value-field] are ignored.""", action='store_true')
    
    parser.add_argument('--no-header', help="""Use this flag if input bedGraph files have no
        header. Standard bedGraph fields are used.""", action='store_true')
    
    parser.add_argument('--add-chr', help="""If this flag is present then 'chr' will be pre-pended
        to sequence names. This is done before any other changes to sequence names, so this
        must be taken into account if giving a dictionary mapping.""", action='store_true')
    
    parser.add_argument('--chr-dict', help="""A dictionary mapping of sequence names found
        in the data to the sequence names, as in "chrom.sizes". The format is as follows:
        '{"key/old":"value/new"}'""", type=json.loads)

    parser.add_argument('--chr-file', help="""A dictionary mapping of sequence names to 
        the sequence names, as in "chrom.sizes", given as a CSV file without header: old,new.
        The file can contain any association, only those in the data will be used.""")

    parser.add_argument('--split-strand', help="""For BED4+ files with strand information
        given by the 5th field as either '+' or '-'. If this option is present, then 
        bedGraph -> bigWig files will be generated for each strand separately. Keywords 
        for forward and reverse strand must be given as '--split-strand FWD REV'.""", 
        nargs=2, type=str, metavar=('FWD', 'REV'))
    
    parser.add_argument('--value-field', help="""For BED4+ files, field name to use as 4th
        field of the bedGraph.""", type=str)

    utils.add_file_options(parser)
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    msg = "[get-bedGraph2bigWig]: {}".format(' '.join(sys.argv))
    logger.info(msg)
    
    # check that stand-alone executable(s) are callable
    programs = ['bedGraphToBigWig']
    utils.check_programs_exist(programs)
    
    # check output path
    if os.path.exists(args.outputDir):
        args.outputDir = os.path.join(args.outputDir, '')
    else:
        msg = "Invalid output path or wrong permission: {}. Terminating.".format(args.outputDir)
        raise OSError(msg)    
    
    # input format - fetch input BED/bedGraph files
    if args.input_format == 'glob': 
        match = r'.*({}|{}.gz).*'.format('|'.join(args.extension, args.extension))
        filenames = list(_glob_re(match, os.listdir(args.inputDir)))
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
            outputFile = os.path.basename(f)
            if outputFile.endswith('gz'):
                outputFile = os.path.splitext(os.path.splitext(outputFile)[0])[0]
            else:
                outputFile = os.path.splitext(outputFile)[0]
            fileMapping[outputFile] = f
    
    # check [--split-strand] keywords
    if agrs.split_strand is not None:
        match = 'fwd|forward|plus|pos|\+|sense|first'
        if re.findall(match, agrs.split_strand[1]):
            msg = """Verify the order of [--split-strand FWD REV]! Keywords were detected
            that might not reflect the correct strand..."""
            logger.warning(msg)
    
    # convert to bigWig directly...
    if args.skip:
        msg = """Using [--skip] and directly converting input files!"""
        for newBed, oldBed in fileMapping.items():
            if not os.path.exists(oldBed):
                msg = "Could not find the bedGraph file: {}. Terminating.".format(oldBed)
                raise FileNotFoundError(msg)
            _convert(oldBed, 
                     os.path.join(args.outputDir, '{}.bw'.format(newBed)),
                     args)
    #... or prepare bedGraph files prior to convert to bigBed.
    else:
        for newBed, oldBed in fileMapping.items():
            if not os.path.exists(oldBed):
                msg = "Could not find the bedGraph file: {}. Terminating.".format(oldBed)
                raise FileNotFoundError(msg)
            if agrs.split_strand is not None:
                for strandStr, strand in zip(agrs.split_strand, ['+', '-']):
                    basename = '{}_{}'.format(newBed, strandStr)
                    filename = _get_bed(oldBed, basename, args, strand=strand)
                    _convert(filename, 
                             os.path.join(args.outputDir, '{}.bb'.format(basename)), 
                             args)
            else:
                filename = _get_bed(oldBed, newBed, args)
                _convert(filename, 
                         os.path.join(args.outputDir, '{}.bb'.format(newBed)), 
                         args)
 

if __name__ == '__main__':
    main()
