#! /usr/bin/env python3

"""Ad hoc script: convert PAS peaks to bigWig fmt
for visualisation.
"""

import sys
import os
import glob
import re
import argparse
import logging
import csv
import json
import pandas as pd

import pbio.misc.shell_utils as shell_utils
import pbio.misc.pandas_utils as pandas_utils
import pbio.misc.logging_utils as logging_utils

import pbio.utils.bed_utils as bed_utils

logger = logging.getLogger(__name__)


BED7_FIELDS = bed_utils.bed6_field_names.copy()
BED7_FIELDS.extend(['count'])
BEDGRAPH_FIELDS = ['seqname', 'start', 'end', 'count']


def get_id(row):
    return '{}:{}-{}'.format(row['seqname'], row['start'], row['end'])


def get_bed_fields(row):
    seqname, coords = row['id'].split(':')
    start, end = coords.split('-')
    row['seqname'] = seqname
    row['start'] = int(start)
    row['end'] = int(end)
    return row


def discard_overlapping_bins(bed, day, rep, strand, args):
    overlaps = bed_utils.get_bed_overlaps(bed, bed)
    problematic = []
    for o in overlaps:
        if o.a_info != o.b_info:
            problematic.append(o.a_info)
            problematic.append(o.b_info)
    problematic = pd.DataFrame(problematic, columns=['id'])
    problematic = problematic.apply(get_bed_fields, axis=1)
    problematic = problematic[['seqname', 'start', 'end', 'id']]
    # save for reference
    output_filename = 'overlap_{}.{}_rep{}.bedGraph'.format(strand, day, rep)
    output_filename = os.path.join(args.tmp, output_filename)
    pandas_utils.write_df(problematic,
                          output_filename,
                          index=False,
                          sep='\t',
                          header=False,
                          do_not_compress=True,
                          quoting=csv.QUOTE_NONE)
    bed = bed.merge(problematic, on=['seqname', 'start', 'end'], how='left')
    bed = bed.loc[~bed.index.isin(bed.dropna(subset=['id_y']).index)]
    return bed


def _convert(bg, bw, args):

    in_files = [bg, args.chrSizes]
    out_files = [bw]
    cmd = "bedGraphToBigWig {} {} {}".format(bg, args.chrSizes, bw)
    
    shell_utils.call_if_not_exists(cmd,
                                   out_files,
                                   in_files=in_files,
                                   overwrite=args.overwrite,
                                   call=True)
    if args.keep:
        return
    try:
        os.remove(bg)
        msg = "Removing: {}".format(bg)
        logger.info(msg)
    except OSError:
        msg = "Could not remove: {}".format(bg)
        logger.info(msg)
        

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Convert PAS peaks to bigWig files by 
        calling the executable program 'bedGraphToBigWig'. PAS peaks files are split into
        forward and reverse, and have the .bed extension.""")

    parser.add_argument('chrSizes', help="The 'chrom.sizes' file for the UCSC database.")

    parser.add_argument('dirloc', help="""The input directory.""")

    parser.add_argument('output', help="""The output directory.""")
    
    parser.add_argument('-cd', '--chr-dict', help="""A dictionary mapping of sequence names found
        in the data to the sequence names, as in "chrom.sizes". The format is as follows:
        '{"key/old":"value/new"}'""", type=json.loads)
    
    parser.add_argument('-cf', '--chr-file', help="""A dictionary mapping of sequence names to 
        the sequence names, as in "chrom.sizes", given as a CSV file without header: old,new.
        The file can contain any association, only those in the data will be used.""")
    
    parser.add_argument('-k', '--keep', help='''If this flag is present, then intermediary bedGraph
        files are not deleted.''', action='store_true')
    
    parser.add_argument('-o', '--overwrite', help='''If this flag is present, then existing files
        will be overwritten.''', action='store_true')

    parser.add_argument('--tmp', help='''tmp directory.''')

    logging_utils.add_logging_options(parser)
    args = parser.parse_args()
    logging_utils.update_logging(args)
    
    # TODO: to discuss with Thiago...
    # library preparation and handling of the peaks
    strand_mapping = {'forward': 'rev',
                      'reverse': 'fwd'}
    ref = '|'.join(strand_mapping.keys())
    
    # check output path
    if os.path.exists(args.output):
        args.output = os.path.join(args.output, '')
    else:
        msg = "Invalid output path or wrong permission: {}. Terminating.".format(args.output)
        raise OSError(msg)
    
    msg = "[peaks_bedGraph2bigWig]: {}".format(' '.join(sys.argv))
    logger.info(msg)
    
    file_count = 0
    for bed in glob.iglob(os.path.join(args.dirloc, '**'), recursive=True):
        if os.path.isfile(bed): 
            strand = re.findall(ref, bed)
            if strand and (os.path.splitext(bed)[-1] == '.bed'):
                # need to get the directory for the rep number, diff day, then fwd and/or rev
                day, rep = os.path.split(bed)[0].split('/')[-1].split('_')
                # read data
                bed_graph = pd.read_csv(bed, 
                                        sep='\t',
                                        header=None,
                                        names=BED7_FIELDS)
                # make sure the first column (seqname) is treated as a string
                bed_graph['seqname'] = bed_graph['seqname'].astype(str)
                
                # there are some overlap on the data...
                # now we ignore these contigs!
                bed_graph['id'] = bed_graph.apply(get_id, axis=1)
                bed_graph = discard_overlapping_bins(bed_graph, day, rep, strand[0], args)
                # to bedGraph
                bed_graph = bed_graph[BEDGRAPH_FIELDS]
                # map chrom names
                if args.chr_dict:
                    for chrom_old, chrom_new in args.chr_dict.items():
                        seqname_m = bed_graph['seqname'] == str(chrom_old)
                        bed_graph.loc[seqname_m, 'seqname'] = str(chrom_new)
                if args.chr_file:
                    chr_map = pd.read_csv(args.chr_file, 
                                          header=None,
                                          index_col=0, 
                                          squeeze=True).to_dict()
                    bed_graph.replace({"seqname": chr_map}, inplace=True)
                    
                # some duplicates...?
                dups = bed_graph.duplicated()
                if dups.any():
                    msg = "Removing duplicate entries from {}".format(bed)
                    logger.warning(msg)
                    bed_graph = bed_graph[~dups]

                # Sort on the chrom field, and then on the chromStart field.
                bed_graph.sort_values(['seqname', 'start'], ascending=[True, True], inplace=True)
                # write bedGraph to disk
                output_filename = 'pas_peaks_{}.{}_rep{}.bedGraph'.format(strand_mapping[strand[0]], day, rep)
                output_filename = os.path.join(args.output, output_filename)
                if os.path.exists(output_filename) and not args.overwrite:
                    msg = "Output file {} already exists. Skipping.".format(output_filename)
                    logger.warning(msg)
                else:
                    file_count += 1
                    pandas_utils.write_df(bed_graph,
                                          output_filename,
                                          index=False,
                                          sep='\t',
                                          header=False,
                                          do_not_compress=True,
                                          quoting=csv.QUOTE_NONE)
            
    msg = "Processed {} files to bedGraph".format(file_count)
    logger.info(msg)
    
    # now converting to bigWig
    for bed in glob.iglob(os.path.join(args.output, '**'), recursive=True):
        if os.path.isfile(bed):
            bigwig = '{}.bw'.format(os.path.splitext(os.path.basename(bed))[0])
            bigwig = os.path.join(args.output, bigwig)
            msg = "Converting {} to {}".format(bed, bigwig)
            logger.info(msg)
            _convert(bed, bigwig, args)


if __name__ == '__main__':
    main()
