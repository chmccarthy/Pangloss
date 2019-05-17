# PanGLOSS


## Requirements
- Python (>2.7.10)
    - BioPython (
- Perl
    - YAML
    - Logger::Simple
    - Parallel::Manager
- R (>3.5.2)
    - Bioconductor (>3.8)
    - ggplot2
    - ggrepel
    - KaryoploteR (>1.8.5)
    - UpSetR
- GeneMark-ES (4.30)        (*)
- TransDecoder (5.0.2)      (*)
- PanOCT (3.23)             (*)
- Exonerate (2.2)
- BLAST+ (2.7.1)
- BUSCO (3.1.0)
- PAML (4.8a)
- MUSCLE (3.8.31)
- InterProScan (5.33-72.0)  (Linux only)
- GOATools

## Dependencies and installation instructions

The installation locations of following programs can be provided to Pangloss in `src/config.ini` or in a custom config file as specified in the example commands above.

- **Exonerate**
   - Exonerate is used to predict gene models based on translated homologs of reference proteins. A continuation of Exonerate is hosted at <a>https://github.com/nathanweeks/exonerate</a>. Exonerate can be installed from source, from apt-get on Linux or through Homebrew on macOS.

- **GeneMark-ES**
  - GeneMark-ES is used for HMM-dependent gene model prediction. macOS and Linux versions of GeneMark-ES are available at <a>http://topaz.gatech.edu/GeneMark/license_download.cgi</a>. GeneMark-ES requires the YAML, Hash::Merge, Logger::Simple and Parallel::ForkManager Perl modules which are all available via cpanm.
  - Pangloss uses a modified version of `get_sequence_from_GTF.pl` provided in `src` rather than the version included with GeneMark-ES.

- **TransDecoder**
  - TransDecoder is used for PVM-dependent gene model prediction. TransDecoder is available from <a>https://github.com/TransDecoder/TransDecoder</a>.

- **yn00**
   - yn00 is required for selection analysis of syntenic ortholog clusters, and is included as part of the PAML software package. PAML source code/binaries are available at <a>http://abacus.gene.ucl.ac.uk/software/paml.html#download</a>.

- **MUSCLE**
   - MUSCLE is required for yn00 analysis of syntenic ortholog clusters. MUSCLE binaries can be found at <a>https://www.drive5.com/muscle/downloads.htm</a>.

- **BUSCO**
    - BUSCO is used to assess gene model set completedness against a set of lineage-specific universal orthologs. Installation instructions for BUSCO are available at <a>https://gitlab.com/ezlab/busco</a>. For completedness analysis of protein sequence data HMMER must also be installed (available from http://hmmer.org/). Note that you need to specify a config.ini file for BUSCO analysis (generally located in $BUSCO_INSTALL_PATH/scripts/../config/) and you need to change the location of HMMsearch to where you have installed the HMMER suite (e.g. /usr/local/bin) in that file.
