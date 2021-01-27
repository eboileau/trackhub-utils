#! /usr/bin/env python3

"""Create hub components.
"""

import sys
import os
import re
import argparse
import logging
import yaml

import trackhub

import phub.utils as utils

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Wrapper to facilitate the creation of
        hub components using the package 'trackhub'.""")

    parser.add_argument('config', help="""The yaml configuration file.""")
    
    parser.add_argument('stagingDir', help="""The local output directory.""")
    #TODO: staging/upload, copy/link source files
    # args for trackhub
    defaultPos = None
    
    utils.add_file_options(parser)
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    msg = "[get-hub]: {}".format(' '.join(sys.argv))
    logger.info(msg)
    
    # check output path
    if os.path.exists(args.stagingDir):
        args.stagingDir = os.path.join(args.stagingDir, '')
    else:
        msg = "Invalid output path or wrong permission: {}. Terminating.".format(args.stagingDir)
        raise OSError(msg)    
    
    # read configuration file
    config = yaml.load(open(args.config), Loader=yaml.FullLoader)
    # first check if we have at least one file format defined
    l = [re.search('^big.+?Files$', k) for k in config.keys()]
    if l.count(None) == len(l):
        msg = """There are no files listed in the configuration file and/or
        the wrong key was used. See the example configuration file provided
        with the package."""
        raise KeyError(msg)
    required_keys = ['hub', 'genomesFile', 'email']
    for fmt in ['Bed', 'Wig']: # supported file format so far...
        files = 'big{}Files'.format(fmt)
        if files in config:
            required_keys.extend([files, 'big{}GlobalSettings'.format(fmt)])
    utils.check_keys_exist(config, required_keys)
    
    # default hub settings    
    short_label = config.get('shortLabel', config['hub'])
    try:
        long_label = config.get('longLabel', 'shortLabel')
    except:
        long_label = config.get('longLabel', config['hub'])

    hub = trackhub.Hub(hub=config['hub'],
                       short_label=short_label,
                       long_label=long_label,
                       email=config['email'],
                       filename='hub.txt')

    genome_kwargs = {}
    if defaultPos:
        genome_kwargs['defaultPos'] = defaultPos
    genome = trackhub.Genome(config['genomesFile'], **genome_kwargs)
    genomesFile = trackhub.GenomesFile(filename='genomes.txt')
    trackDb = trackhub.TrackDb()
    hub.add_genomes_file(genomesFile)
    genomesFile.add_genome(genome)
    genome.add_trackdb(trackDb)

    # standard settings - individual tracks and/or tracks grouped with superTrack container
    # check if we have superTracks, else add tracks
    superTracks = {}
    if 'superTracks' in config.keys():
        for superTrack, tracks in config['superTracks'].items():
            superTracks[superTrack] = tracks
            defaultSuperTrackSettings = {'shortLabel': superTrack, 'longLabel': superTrack}
            settings = config['superTrackSettings'].get(superTrack, defaultSuperTrackSettings)
            short = settings.get("shortLabel", superTrack)
            long = settings.get("longLabel", superTrack)
            s = trackhub.SuperTrack(
                name=superTrack,
                short_label=short,
                long_label=long,
            )    
            exec('super{}=s'.format(superTrack))
            exec('trackDb.add_tracks(super{})'.format(superTrack))
        
    # first add tracks associated with superTracks, then add remaining tracks if any
    for superTrack, tracks in superTracks.items():
        for track in tracks:
            for fileType, globalSettings, fileSettings in zip(['bigBedFiles', 'bigWigFiles'], 
                                                            ['bigBedGlobalSettings', 'bigWigGlobalSettings'],
                                                            ['bigBedFileSettings', 'bigWigFileSettings']):
                if track in config[fileType]:
                    source = config[fileType][track]
                    gSettings = config[globalSettings].copy()
                    trackType = gSettings.pop('trackType')
                    fSettings = config[fileSettings].get(track, None)
                    if fSettings is not None:
                        finalSettings = {**gSettings, **fSettings}
                    else:
                        finalSettings = gSettings
                    # define track
                    track = trackhub.Track(
                        name=track,
                        source=os.path.basename(source), # track
                        url=os.path.basename(source),
                        tracktype=trackType)
                    for setting, value in finalSettings.items():
                        exec("track.add_params({}='{}')".format(setting, value))
                    # add to superTrack
                    exec('super{}.add_tracks(track)'.format(superTrack))

    # now process tracks not attached to a superTrack      
    processed = sum(superTracks.values(), [])  
    for fileType, globalSettings, fileSettings in zip(['bigBedFiles', 'bigWigFiles'], 
                                                    ['bigBedGlobalSettings', 'bigWigGlobalSettings'],
                                                    ['bigBedFileSettings', 'bigWigFileSettings']):

        for track in config[fileType]:
            if track not in processed:
                source = config[fileType][track]
                gSettings = config[globalSettings].copy()
                trackType = gSettings.pop('trackType')
                fSettings = config[fileSettings].get(track, None)
                if fSettings is not None:
                    finalSettings = {**gSettings, **fSettings}
                else:
                    finalSettings = gSettings
                # define track
                track = trackhub.Track(
                    name=track,
                    source=os.path.basename(source), # track
                    url=os.path.basename(source),
                    tracktype=trackType)
                for setting, value in finalSettings.items():
                    exec("track.add_params({}='{}')".format(setting, value))
                # add
                trackDb.add_tracks(track)    
    
    # TODO: files need to excist to be linked to the staging dir, currently only
    # the hub files/structure is created, so this ends with ValueError: target {} not found.
    
    # stage the hub (no upload)
    trackhub.upload.stage_hub(hub, staging=args.stagingDir)

        
if __name__ == '__main__':
    main()
    
