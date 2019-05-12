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

## Installion instructions

- PAML
   - PAML source code/binaries are available at <a>http://abacus.gene.ucl.ac.uk/software/paml.html#download</a>.

- MUSCLE
   - MUSCLE binaries can be found at <a>https://www.drive5.com/muscle/downloads.htm</a>. MUSCLE is required for yn00 analysis of translated protein families.

- BUSCO
    - Installation instructions for BUSCO are available at <a>https://gitlab.com/ezlab/busco</a>. For completedness analysis of protein sequence data HMMER must also be installed (available from http://hmmer.org/). Note that you need to specify a config.ini file for BUSCO analysis (generally located in $BUSCO_INSTALL_PATH/scripts/../config/) and you need to change the location of HMMsearch to where you have installed the HMMER suite (e.g. /usr/local/bin) in that file.
