# Description

This repository contains scripts, wrappers or snippets of code useful for converting
format ( bigBed, bigWig, *etc.* ) and creating track hubs.

## Getting Started

First copy the application binaries

```
cd scripts 
./copy_ucsc_exe
```

and grab chrom.sizes from UCSC using the chosen assembly *e.g.*

```
mkdir ../local
./fetchChromSizes hg38 > ../local/hg38.chrom.sizes
```

Perform any operation necessary to prepare the files that will be used to create tracks, such as

```
# set uniform color in BED files...
./set_bed_color BED1.bed BED2.bed.gz 265,165,0
# ... remove selected contigs, adjust score, etc.
header=$(zcat BED1.bed.gz | sed -n '1p') && zcat BED1.bed.gz | tail $1 -n+2 | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$3,$4,100,$6,$7,$8,$9,$10,$11,$12;}' | awk '$1!~"AEMK|FPK"' | sed "1s/^/$header\n/g" | gzip > BED2.bed.gz
# reheader BAM files to match UCSC contig format...
./reheader_bam DIR1 DIR2 hg38
# ... or run any other custom script to convert GTF, BED6 to BED12, or bedGraph.
```

Prepare the indexed binary files (bigBed, bigWig) by running one of the following `get_bed2bigBed.py`, `get_bam2bigWig.py`, or `get_bedGraph2bigWig.py`.
These scripts are wrapper to Genome Browser application binaries or `deepTools` executables such as `bamCoverage`. For more information, use the 
`--help\-h` option.

TODO:
- clean `get_bed2bigBed.py`, `get_bam2bigWig.py`, and `get_bedGraph2bigWig.py`.
- go back to pyproc and rewrite `create_bigBed_tracks.py` such that it is only for rpbp output, which can be then used with `get_bed2bigBed.py` here. This
splits the code in 2, but then each part is contained, 1 for rpbp, 2 for all.
- for `get_bed2bigBed.py`, we do not use pbio, and/or assume bed fields, etc.
- then , create base python script to create trackDB (see bash script) as before (using args?), then add additional option from config
- can we call external script in python such that rackhub classes are seen??????

- option: add gtf2bed12, bed62bed12, etc. or pep tobed, etc....???!



### Prerequisites

Pinned version of selected dependencies are listed in the `requirements.txt` file for reproducible installation,
in particular the Python package `trackhub` and `deepTools`.

A number of UCSC Genome Browser application binaries for stand-alone use are also required and can be copied
using the script `copy_ucsc_exe`.


### Installation

To install the local VCS project in development mode, use the `--editable` or `-e` option, otherwise
this flag can be ignored. 

To install `phub` and dependencies, first create a virtual environment:
 
```
python3 -m venv /path/to/virtual/environment
```

For information about Python virtual environments, see the [venv](https://docs.python.org/3/library/venv.html) documentation.
To activate the new virtual environment and install `phub`:

```
# Activate the new virtual environment.
source /path/to/virtual/environment/bin/activate

# If necessary, upgrade pip and wheel or additional packages (such as setuptools if installing in editable mode).
pip install --upgrade pip setuptools wheel

# Clone the git repository
git clone https://github.com/eboileau/trackhub-utils
cd trackhub-utils

# The period is required, it is the local project path (trackhub-utils)
pip --verbose install -r requirements.txt [-e] . 2>&1 | tee install.log

```

#### Anaconda installation

The package can also be installed within an [anaconda](https://www.continuum.io/) environment. 

```
# Create the anaconda environment.
conda create -n my_new_environment python=3.6 anaconda

# Activate the new environment.
source activate my_new_environment

# Clone the git repository
git clone https://github.com/eboileau/trackhub-utils
cd trackhub-utils

pip --verbose install -r requirements.txt [-e] . 2>&1 | tee install.log
```

## Uninstallation

To remove the `phub` package:

```
pip uninstall phub
```

If the package is installed in a dedicated virtual environment, this environment can also be cleared or removed.

## Running the tests

## Contributing

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

Some scripts/utils are adapted from the `pybio-utils` package, authored by Brandon Malone, and maintained by Etienne Boileau.


