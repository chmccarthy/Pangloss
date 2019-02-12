# -*- coding: utf-8 -*-
"""
PanGLOSS: A pipeline for pangenome analysis of microbial eukaryotes.

Requirements:
    - Python (written for 2.7.x)
        - BioPython (1.73)
    - Exonerate (>2.2)
    - GeneMark-ES (>4.30)
    - TransDecoder (>5.0.2)
    - BLAST+ (>2.7.1)
    -

Recent changes:
    v0.2.0 (February 2019)


    v0.1.0 (January 2019)
    - Added config file and command line parsers.
    - Rewrote PanGuess and how it's handled from master script.
    - Created master script based on old pangenome pipelines.



Written by Charley McCarthy, Genome Evolution Lab, Department of Biology,
Maynooth University in 2017-2019 (Charley.McCarthy@nuim.ie).
"""

import logging
import sys

from argparse import ArgumentParser
from ConfigParser import SafeConfigParser
from datetime import datetime
from glob import glob

from PanGLOSS import BLASTAll, PanGuess, PanOCT, QualityCheck


def PanGuessHandler(genomelist, workdir, ref, exon_cov, gm_branch, td_potenial, td_len, cores=None, skip=False):
    """
    Runs PanGuess from master script.
    
    Arguments taken from Gene_model_prediction section of config file as follows:
        genomelist   = List of strain genomes specified by genomes_list.
        workdir      = Working directory for prediction given by prediction_dir.
        ref          = FASTA file of reference protein set given by reference_proteins.
        exon_cov     = Exonerate sequence coverage cutoff given by
                       exonerate_sequence_coverage (int).
        gm_branch    = Option for fungal-specific branch-site prediction model for
                       GeneMark-ES given by genemark_fungal_model (0 or 1).
        td_potential = Cutoff for TransDecoder coding potenial score given by
                       trans_coding_potential (int).
        td_len       = Amino acid sequence length cutoff for TransDecoder given by
                       trans_aa_length (int).

    Arguments activated by command line flags:
        skip         = Skip Exonerate gene model predictions.

    Find some way to incorporate exon_cov and td_potential into final product!
    """
    # Generate list of genomes from user-provided genome list file.
    logging.info("Master: Parsing genome list.")
    genomes = [line.strip("\n") for line in open(genomelist)]
    
    # Create working directory if not present.
    logging.info("Master: Building working directory for gene model prediction.")
    PanGuess.MakeWorkingDir(workdir)

    #
    if not skip:
        logging.info("Master: Building working directory for gene model prediction.")
        PanGuess.BuildRefSet(workdir, ref)
    
    # Loop over each genome and carry out gene model prediction.
    for genome in genomes:
        # Make tag from genome name.
        tag = genome.split(".")[0].split("/")[1]
        tags.append(tag)
        logging.info("Master: Running gene model prediction for {0}.".format(tag))

        if not skip:
            # Run prediction using Exonerate.
            cmds = PanGuess.BuildExonerateCmds(workdir, genome)
            exonerate_genes = PanGuess.RunExonerate(cmds, cores)
        
            # Order gene models predicted via Exonerate by Contig ID: Location.
            logging.info("Master: Sorting gene model predictions by genomic location.")
            exonerate_genes.sort(key=lambda x: (x.contig_id, x.locs[0]))
        
            # Extract genomic attributes from Exonerate gene model set.
            exonerate_attributes = PanGuess.GetExonerateAttributes(exonerate_genes, tag)

        else:
            logging.info("Master: Skipping gene model prediction via Exonerate (--no_exonerate enabled).")
            exonerate_genes = None
        
        # Run prediction using GeneMark-ES.
        logging.info("Master: Running gene model prediction for {0} using GeneMark-ES.".format(genome))
        genemark_gtf = PanGuess.RunGeneMark(genome, gm_branch, cores)
        
        # Convert GeneMark-ES GTF file into a more PanOCT-compatible version.
        logging.info("Master: Converting GeneMark GTF data to attribute data.")
        genemark_attributes = PanGuess.GeneMarkGTFConverter(genemark_gtf, tag)
        
        # Merge unique gene model calls between Exonerate and GeneMark-ES.
        if not skip:
            logging.info("Master: Merging Exonerate and GeneMark-ES gene calls.")
            merged_attributes = PanGuess.MergeAttributes(exonerate_attributes, genemark_attributes)
        else:
            merged_attributes = genemark_attributes
            del genemark_attributes

        # Clean up GeneMark-ES files and folders.
        logging.info("Master: Tidying up GeneMark-ES temporary files.")
        PanGuess.MoveGeneMarkFiles(workdir, genome)
        
        # Extract NCRs into list.
        logging.info("Master: Extracting non-coding regions of {0} for TransDecoder analysis.".format(genome))
        noncoding = PanGuess.ExtractNCR(merged_attributes, genome)
        
        # Run TransDecoder on NCRs.
        logging.info("Master: Running TransDecoder on non-coding regions of {0}.".format(genome))
        tdir = PanGuess.RunTransDecoder(noncoding, workdir, genome, td_len)
        
        # Move TransDecoder files.
        logging.info("Master: Tidying up TransDecoder temporary files.")
        PanGuess.MoveTransDecoderFiles(tdir)
        
        # Extract TransDecoder attributes.
        logging.info("Master: Converting TransDecoder GTF data to attribute data.")
        trans_attributes = PanGuess.TransDecoderGTFToAttributes(tdir, tag)
        
        # Merge TransDecoder calls into the Exonerate + GeneMark-ES set.
        logging.info("Master: Merging all remmaining gene calls for {0}.".format(genome))
        full_attributes = PanGuess.MergeAttributes(merged_attributes, trans_attributes)
        
        # Write out gene set, protein set and attributes set.
        logging.info("Master: Writing out datasets for {0}.".format(genome))
        PanGuess.ConstructGeneModelSets(full_attributes, exonerate_genes, workdir, genome, tag)
        
        # Compress temporary folders and finish up.
        logging.info("Master: Compressing temporary folders for {0}.".format(genome))
        PanGuess.TarballGenePredictionDirs(workdir, genome)
        logging.info("Master: Finished gene model predction for {0}.".format(genome))


def QualityCheckHandler(tags, queries, cores=None):
    """
    Search a user-provided set of genes of dubious-quality (i.e. pseudogenes, transposable elements or
    transposons &c.) against predicted gene model sets and filter out sufficiently similar genes in the latter.

    Arguments:
        gene_sets   = List of strains in analysis (easy access to all files associated with a strain).
        queries     = Set of genes (protein sequences, in fact) to search against all gene model sets.
    """
    QualityCheck.BuildMakeBLASTDBs(tags, cores)
    blasts = QualityCheck.QCBLAST(queries, tags, cores)
    QualityCheck.RemoveDubiousCalls(blasts, tags)


def BUSCOHandler():
    pass


def BLASTAllHandler(tags, evalue=0.0001, cores=None):
    """
    Runs all-vs.-all BLASTp search of gene model dataset as required for PanOCT.
    """
    BLASTAll.ConcatenateDatasets(tags)
    blasts = BLASTAll.BLASTAll(evalue, cores)
    BLASTAll.MergeBLASTsAndWrite(blasts)


def PanOCTHandler(fasta_db, attributes, blast, tags, **kwargs):
    """
    Runs PanOCT.
    """
    PanOCT.RunPanOCT(fasta_db, attributes, blast, tags, **kwargs)
    PanOCT.PanOCTOutputHandler()


def IPSHandler():
    pass


### Parser functions. ###

def CmdLineParser():
    """
    Create and return a configuration file parser.
    """
    # Create our argument parser
    ap = ArgumentParser(description="Pan-genome analysis of microbial eukaryotes.")

    # The argument group for the gene model prediction step.
    pred = ap.add_mutually_exclusive_group(required=True)

    # Add arguments for gene model prediction.
    pred.add_argument("--pred", action="store_true", help="Perform gene model prediction.")
    pred.add_argument("--pred_only", action="store_true", help="Only perform gene model prediction. "
                      "Allows users to generate their own all-vs.-all BLASTp data (e.g. on HPC environments), "
                      "and rerun full analysis later on.")
    pred.add_argument("--no_pred", action="store_true", help="Skip gene model prediction.")

    # Add argument to skip Exonerate steps during gene model prediction.
    ap.add_argument("--no_exonerate", action="store_true", help="Skip gene model prediction via Exonerate.")

    # Add arguments for optional quality control analyses.
    ap.add_argument("--qc", action="store_true", help="Perform quality check on predicted gene model sets.")
    ap.add_argument("--busco", action="store_true", help="Perform BUSCO analysis on predicted gene model sets.")

    # Add argument for skipping all-vs.-all BLASTp step (usually faster to generate data elsewhere).
    ap.add_argument("--no_blast", action="store_true", help="Skip all-vs.-all BLASTp step for PanOCT.")

    # Add mandatory positional argument for path to config file.
    ap.add_argument("CONFIG_FILE", help="Path to PanGLOSS configuration file.")

    # Parse and return args.
    args = ap.parse_args()
    return args


def ConfigFileParser():
    """
    Create and return a configuration file parser.
    """
    cp = SafeConfigParser()
    return cp


### Main function. ###

def main():
    """
    """
    # Create logfile and assign it to all child modules.
    start_time = datetime.now()
    logging.basicConfig(filename="PanGLOSS_Run_{0}.log".format(str(start_time).replace(" ", "_")),
                        level=logging.INFO, format="%(asctime)s: %(levelname)s: %(message)s")

    # Create parser objects.
    ap = CmdLineParser()
    cp = ConfigFileParser()

    # Read config file.
    logging.info("Master: Read config file.")
    cp.read("config.txt")

    # Unless disabled, parse arguments for PanGuess and run gene model prediction.
    if ap.pred:
        panguess_args = []
        logging.info("Master: Performing gene prediction steps using PanGuess.")
        for arg in cp.items("Gene_model_prediction"):
            panguess_args.append(arg[1])
        if ap.no_exonerate:
            panguess_args.append(True)
        PanGuessHandler(*panguess_args)
        logging.info("Master: Gene prediction finished.")
    else:
        logging.info("Master: Skipped gene prediction steps (--nopred enabled).")

    # If enabled, check gene sets against user-provided sets of dubious genes, or transposable elements, &c.
    if ap.qc:
        logging.info("Master: Performing gene model QC using QualityCheck.")
        qc_args = [[i.split(".")[0] for i in glob("*.faa")]]
        for arg in cp.items("Quality_check"):
            if arg[1]:
                qc_args.append(arg[1])
        QualityCheckHandler(*qc_args)
        logging.info("Master: QC analysis finished.")
    else:
        logging.info("Master: Skipped gene model QC (--qc not enabled).")

    # If enabled, run BUSCO completedness analysis of all gene model sets.
    if ap.busco:
        pass

    # Allow program to finish after gene prediction and (optionally) QC/BUSCO if --pred_only is enabled.
    if ap.pred_only:
        logging.info("Master: Finishing PanGLOSS (--pred_only enabled). To run remaining steps with your own "
                     "all-vs.-all BLAST data, run PanGLOSS with the --no_pred and --no_blast flags.")
        logging.info("Master: PanGLOSS finished in {0}.".format(str(datetime.now() - start_time)))
        sys.exit(0)

    # Run all-vs.-all BLASTp, unless --no_blast is enabled (i.e., user provides own blast file).
    if not ap.no_blast:
        logging.info("Master: Performing all-vs.-all BLASTp searches for entire dataset.")
        blast_args = []
        for arg in cp.items("BLASTAll_settings"):
            if arg[1]:
                blast_args.append(arg[1])
        BLASTAllHandler(*blast_args)
        logging.info("Master: All-vs.-all analysis finished.")
    else:
        logging.info("Master: Skipping all-vs.-all BLASTp searches (--no_blast enabled).")

    # Run PanOCT on full dataset.
    panoct_default_args = []
    panoct_extra_args = []
    for arg in cp.items("PanOCT_settings"):
        if arg[1]:
            panoct_default_args.append(arg[1])

    if panoct_extra_args:
        pass
    else:
        PanOCTHandler(panoct_default_args)


if __name__ == "__main__":
    main()

