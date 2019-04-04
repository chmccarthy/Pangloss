# -*- coding: utf-8 -*-
"""
PanGLOSS: A pipeline for pan-genome analysis of microbial eukaryotes.

Dependencies (version tested) (* = required):
    - Python (2.7.10)           (*)
        - BioPython (1.73)      (*)
    - Perl                      (*)
        - YAML                  (for GeneMark-ES)
        - Logger::Simple        (for GeneMark-ES)
        - Parallel::Manager     (for GeneMark-ES)
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
    - R (>3.5.2)
        - Bioconductor (>3.8)
        - KaryoploteR (>1.8.5)

To-do:
    - BUSCO assessment of gene set completedness.
    - Annotation of pan-genomes using InterProScan.
    - Assess ancestry of pan-genome complements using BLASTp.
    - Improve logging.
    - Improve config file.
        - Add check for dependencies based on paths in file.
        - Add in section for PanOCT parameters.
        - Add in section for yn00 parameters.


Recent changes:
    v0.6.0 (March-April 2019)
    -
    -


    v0.5.0 (February 2019)
    - Made fixes in PanGuess.
    - Added karyotype plotting of PanOCT clusters for all genomes using KaryoploteR (see Karyotype.py and .R files).

    v0.4.0 (February 2019)
    - Added in yn00 module for running selection analysis on nucleotide sequences.
        - A bug in Biopython prevents yn00 from running if a dot is in the sequence name,
          will incorporate a workaround into future versions until bug is fixed.
    - Added functions in Tools for yn00 analysis workflow.
    - Added --yn00 flag to command line.

    v0.3.0 (February 2019)
    - PanGLOSS now requires Biopython >1.73 for correct handling and parsing of SearchIO objects in BLASTAll.
    - Added in BLASTAll module to handle (optional) all-vs.-all BLASTp searches for PanOCT.
    - Added in PanOCT module to handle PanOCT analysis with default parameters.
    - Redesigned workflow to allow user to solely carry out gene model prediction.
    - Redesigned workflow to allow user to skip prediction and BLASTAll steps if specified.
    - Redesigned main function to reflect optional arguments/workflows.
    - Improved logging.

    v0.2.0 (February 2019)
    - Added in QualityCheck module to handle (optional) filtering for pseudogenes.
    - Incorporated logging.

    v0.1.0 (January 2019)
    - Constructed basic version of PanGLOSS-compatible config file.
    - Added config file and command line parsers.
    - Rewrote PanGuess as module, and changed how it's handled from master script.
    - Created master script based on old pangenome pipelines from 2017-18.

Written by Charley McCarthy, Genome Evolution Lab, Department of Biology,
Maynooth University in 2017-2019 (Charley.McCarthy@nuim.ie).
"""

import logging
import sys
from argparse import ArgumentParser
from ConfigParser import SafeConfigParser
from datetime import datetime
from glob import glob
from PanGLOSS import BLASTAll, GO, PAML, PanGuess, PanOCT, QualityCheck, Karyotype


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

    # Build reference protein directory for Exonerate (unless --no_exonerate is enabled).
    if not skip:
        logging.info("Master: Building working directory for gene model prediction.")
        PanGuess.BuildRefSet(workdir, ref)
    
    # Loop over each genome and carry out gene model prediction.
    for genome in genomes:
        # Make tag from genome name.
        tag = genome.split(".")[0].split("/")[1]
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
            pass
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
    # Build BLAST DB, run QC searches against DB and filter out any dubious gene calls.
    logging.info("Master: Running QualityCheckHandler.")
    QualityCheck.BuildMakeBLASTDBs(tags, cores)
    blasts = QualityCheck.QCBLAST(queries, tags, cores)
    QualityCheck.RemoveDubiousCalls(blasts, tags)


def BUSCOHandler():
    """
    Run BUSCO on all gene sets in a pan-genome dataset. Returns a text file detailing completedness
    """
    pass


def BLASTAllHandler(tags, evalue=0.0001, cores=None):
    """
    Runs all-vs.-all BLASTp search of gene model dataset as required for PanOCT.
    """
    # Concatenate all protein sequence datasets together, BLAST them against themselves,
    # pool all farmed results together and write output (in tabular format) to file.
    logging.info("Master: Running BLASTAllHandler.")
    BLASTAll.ConcatenateDatasets(tags)
    blasts = BLASTAll.BLASTAll(evalue, cores)
    BLASTAll.MergeBLASTsAndWrite(blasts)


def PanOCTHandler(fasta_db, attributes, blast, tags, gaps=False, **kwargs):
    """
    Runs PanOCT and does some post-run cleanup and sequence extraction.
    """
    # Run PanOCT with provided files (and optional additional arguments.
    logging.info("Master: Running PanOCTHandler.")
    #PanOCT.RunPanOCT(fasta_db, attributes, blast, tags, **kwargs)

    # If enabled, try to fill potential gaps in syntenic clusters within pangenome using BLAST+ data.
    if gaps:
        logging.info("Master: Running gap filling method.")
        PanOCT.FillGaps(blast, "matchtable.txt", fasta_db, tags)

    # Move all output from PanOCT into dedicated subfolder, and extract syntenic clusters to their own subfolder.
    #PanOCT.PanOCTOutputHandler()
    #PanOCT.GenerateClusterFASTAs()


def GOHandler():
    """
    Run GO-slim enrichment analysis on pangenome datasets using GOATools.
    """
    # Make GO folder, unless one already exists.
    GO.MakeWorkingDirs()

    # Generate dictionary for IPS annotation data.
    annos = GO.GenerateAnnoDict("yarr_test.tsv")

    # Generate GO associations and populations files.
    GO.GenerateAssociations(annos)
    GO.GeneratePopulations(annos, "matchtable.txt")

    # Run map_to_slim.py from GOATools.
    GO.GenerateSlimData("go/associations.txt", "go.obo", "goslim_generic.obo")

    # Run enrichment analysis of core and accessory genomes against full pangenome dataset.
    GO.CoreEnrichment("go.obo", "go/core_pop.txt", "go/full_pop.txt", "go/pangenome_slim.txt")
    GO.AccessoryEnrichment("go.obo", "go/acc_pop.txt", "go/full_pop.txt", "go/pangenome_slim.txt")


def PAMLHandler():
    """
    Run Yn00 on core and accessory gene model clusters.
    """
    seqs = PAML.TranslateCDS()
    alignment = PAML.MUSCLEAlign(seqs)
    PAML.PutGaps(alignment, seqs)
    PAML.RunYn00(seqs)


def KaryoploteRHandler():
    """
    Generates chromosomal plots of core and accessory gene models for each genome in a dataset, similar to
    the Ruby program PhenoGram but with way less overhead.
    """
    # Make lengths and karyotypes files (don't think these can be passed as objects to R without a lot of effort).
    Karyotype.GenerateContigLengths("genomes")
    if ap.gaps:
        Karyotype.GenerateKaryotypeFiles("allatt.db", "panoct/refined_matchtable.txt")
    else:
        Karyotype.GenerateKaryotypeFiles("allatt.db", "panoct/matchtable.txt")

    # Pass required files to KaryPloteR and run R script.
    Karyotype.KaryoPloteR("./panoct_tags.txt", "./karyotypes.txt", "./genomes/lengths.txt")


### Parser functions. ###

def CmdLineParser():
    """
    Create and return a configuration file parser.
    """
    # Create our argument parser.
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

    # Add argument for skipping all-vs.-all BLASTp step (usually faster to generate data elsewhere).
    ap.add_argument("--no_blast", action="store_true", help="Skip all-vs.-all BLASTp step for PanOCT.")

    # Add arguments for optional quality control analyses.
    ap.add_argument("--qc", action="store_true", help="Perform quality check on predicted gene model sets.")
    ap.add_argument("--busco", action="store_true", help="Perform BUSCO analysis on predicted gene model sets.")

    # Add argument for gap filling in PanOCT-dervied pangenome.
    ap.add_argument("--fillgaps", action="store_true", help="Attempt to fill potential gaps in syntenic clusters.")

    # Add arguments for annotation and GO-enrichment analysis.
    ap.add_argument("--ips", action="store_true", help="Perform InterProScan analysis of gene model sets. NOTE:"
                                                       " Do not enable this option on non-Linux operating systems,"
                                                       " InterProScan is not supported on these systems.")
    ap.add_argument("--go", action="store_true", help="Perform GO-slim enrichment analysis using GOATools.")

    # Add argument for selection analysis using yn00.
    ap.add_argument("--yn00", action="store_true", help="Perform selection analysis on core and accessory gene "
                                                        "families using yn00.")
    # Add argument to produce karyotype plots.
    ap.add_argument("--karyo", action="store_true", help="Generate karyotype plots for all genomes in a "
                                                         "database based on PanOCT results.")

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
    Main function.
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
    if any([ap.pred, ap.pred_only]):
        panguess_args = []
        logging.info("Master: Performing gene prediction steps using PanGuess.")
        for arg in cp.items("Gene_model_prediction"):
            panguess_args.append(arg[1])
        if ap.no_exonerate:
            panguess_args.append(True)
        PanGuessHandler(*panguess_args)
        logging.info("Master: Gene prediction finished.")

        # Allow program to finish after gene prediction and (optionally) QC/BUSCO if --pred_only is enabled.
        if ap.pred_only:
            logging.info("Master: Finishing PanGLOSS (--pred_only enabled). To run remaining steps with your own "
                         "all-vs.-all BLAST data, run PanGLOSS with the --no_pred and --no_blast flags.")
            logging.info("Master: PanGLOSS finished in {0}.".format(str(datetime.now() - start_time)))
            sys.exit(0)
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
    if ap.fillgaps:
        panoct_default_args.append(True)
    if panoct_extra_args:
        pass
    else:
        PanOCTHandler(*panoct_default_args)

    # If enabled, run InterProScan analysis on entire dataset.
    if ap.ips:
        if sys.platform != "linux":
            print "InterProScan is not supported on non-Linux operating systems. Cannot run InterProScan analysis."
            print "See https://github.com/ebi-pf-team/interproscan/wiki for more information."
            sys.exit(0)
        else:
            pass

    # If enabled, run GO-slim enrichment analysis on core and accessory datasets using GOATools.
    if ap.go:
        GOHandler()
        pass

    # If enabled, run selection analysis using yn00.
    #if ap.yn00:
    #    logging.info("Master: Performing selection analysis using yn00.")
    #    PAMLHandler()

    # If enabled, generate karyotype plots for all strain genomes in pangenome dataset.
    if ap.karyo:
        logging.info("Master: Generating karyotype plots for all genomes in dataset.")
        KaryoploteRHandler()


if __name__ == "__main__":
    main()

