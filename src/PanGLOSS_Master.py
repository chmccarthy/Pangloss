# -*- coding: utf-8 -*-
"""
Pangloss: A pipeline for pan-genome analysis of microbial eukaryotes.

Written by Charley McCarthy, Genome Evolution Lab, Department of Biology,
Maynooth University between 2017-2019 (Charley.McCarthy@nuim.ie).

See Pangloss/README.md for more information.

To-do:
    - Improve logging.
    - Seriously improve config file!
    - Test IPS pipeline fully.

Recent changes:
    v0.8.0 (May 2019)
    - Rejigged when exactly amino acid, nucleotide and attribute datasets are concatenated.
    - Loads of changes to R scripts.
    - Fully implemented pangenome refinement using "gap finding" method.

    v0.7.0 (May 2019)
    - Changed --fillgaps to --refine to reflect output matchtable.
    - Added percentage score threshold of ≥90% to Exonerate gene model prediction in place of sequence coverage.
    - Added BUSCO assessment of gene model set completedness.
    - Added option for running InterProScan analysis of dataset within PanGLOSS.
    - Fully implemented yn00 selection analysis and summary generation for full datasets.
    - Moved cores check out of modules and into master script.
    - Fully incorporated paths from config file into script and set default config file.

    v0.6.0 (May 2019)
    - Added in command flag for disabling PanOCT (mostly for debugging purposes).
    - Added method for refining pangenome construction based on BLAST+ data used by PanOCT.
    - Added in bar chart plotting of pangenome cluster sizes using ggplot and a function for estimating true
      pangenome size using Chao's lower bound estimation method.
    - Added in ring chart plotting of pangenome component sizes using ggplot2 and ggrepel (see Size.py and RingChart.R
      files).
    - Added in UpSet plotting of distribution of syntenic orthologs within accessory genome using UpSetR (see UpSet.py
      and .R files).
    - Slight change to what goes into karyotype input file, now includes gene names. Done this with a view to a making
      a form of gene order/synteny chart in a possible future version.
    - Changed how subdirectories are created everywhere.

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
"""

import logging
import os
import sys
from ConfigParser import SafeConfigParser
from argparse import ArgumentParser
from datetime import datetime
from glob import glob

from PanGLOSS import BLASTAll, BUSCO, GO, Karyotype, PAML, PanGuess, PanOCT, QualityCheck, Size, UpSet
from PanGLOSS.Tools import ConcatenateDatasets


def PanGuessHandler(ex_path, gm_path, tp_path, tl_path,
                    genomelist, workdir, ref, gm_branch, td_len, cores=None, skip=False):
    """
    Runs PanGuess from master script.

    Paths taken from config file:
        ex_path      = Exonerate path.
        gm_path      = GeneMark-ES path.
        tp_path      = TransDecoder.Predict path.
        tl_path      = TransDecoder.LongOrfs path.
    
    Arguments taken from Gene_model_prediction section of config file as follows:
        genomelist   = List of strain genomes specified by genomes_list.
        workdir      = Working directory for prediction given by prediction_dir.
        ref          = FASTA file of reference protein set given by reference_proteins.
        gm_branch    = Option for fungal-specific branch-site prediction model for
                       GeneMark-ES given by genemark_fungal_model (0 or 1).
        td_len       = Amino acid sequence length cutoff for TransDecoder given by
                       trans_aa_length (int).
        cores        = Number of threads to run predictions on.

    Arguments activated by command line flags:
        skip         = Skip Exonerate gene model predictions.
    """
    # If user doesn't specify cores in config file, just leave them with one free.
    if not cores:
        cores = str(mp.cpu_count() - 1)

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
            cmds = PanGuess.BuildExonerateCmds(workdir, ex_path, genome)
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
        genemark_gtf = PanGuess.RunGeneMark(genome, gm_path, gm_branch, cores)
        
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
        tdir = PanGuess.RunTransDecoder(noncoding, tp_path, tl_path, workdir, genome, td_len)
        
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

    ConcatenateDatasets(genomelist)


def QualityCheckHandler(sets, queries, cores=None):
    """
    Search a user-provided set of genes of dubious-quality (i.e. pseudogenes, transposable elements or
    transposons &c.) against predicted gene model sets and filter out sufficiently similar genes in the latter.

    Arguments:
        gene_sets   = List of strains in analysis (easy access to all files associated with a strain).
        queries     = Set of genes (protein sequences, in fact) to search against all gene model sets.
    """
    # Build BLAST DB, run QC searches against DB and filter out any dubious gene calls.
    logging.info("Master: Running QualityCheckHandler.")
    QualityCheck.BuildMakeBLASTDBs(sets, cores)
    blasts = QualityCheck.QCBLAST(queries, sets, cores)
    QualityCheck.RemoveDubiousCalls(blasts, sets)


def BUSCOHandler(buscopath, lineagepath, tags):
    """
    Run BUSCO on all gene sets in a pangenome dataset. Returns a text file detailing completedness of each
    gene model set in pangenome dataset.
    """
    BUSCO.RunBUSCO(buscopath, lineagepath, tags)
    pass


def BLASTAllHandler(tags, cores=None):
    """
    Runs all-vs.-all BLASTp search of gene model dataset as required for PanOCT. BLASTp searches "parallelized"
    via subprocessing and cStringIO magic. Can be skipped from command-line, and in general it might be better
    for the user to run all-vs.-all BLASTp searches on some kind of HPC server if possible.

    Arguments:
        tags   = List of strains in analysis (easy access to all files associated with a strain).
        evalue = E-value cutoff for BLASTp searches (default is 10^-4).
        cores  = Number of BLASTp searches to run simulatenously (default will be available cores - 1).
    """
    # Concatenate all protein sequence datasets together, BLAST them against themselves,
    # pool all farmed results together and write output (in tabular format) to file.
    logging.info("Master: Running BLASTAllHandler.")
    ConcatenateDatasets(tags)
    blasts = BLASTAll.BLASTAll(cores)
    BLASTAll.MergeBLASTsAndWrite(blasts)


def PanOCTHandler(fasta_db, attributes, blast, tags, gaps=False, **kwargs):
    """
    Runs PanOCT, refines initial construction if enabled and does some post-run cleanup and sequence extraction.

    Arguments:
        gene_sets   = List of strains in analysis (easy access to all files associated with a strain).
        queries     = Set of genes (protein sequences, in fact) to search against all gene model sets.
    """
    # Run PanOCT with provided files (and optional additional arguments.
    logging.info("Master: Running PanOCTHandler.")
    if not os.path.isfile(fasta_db):
        ConcatenateDatasets("genomes/genomes.txt")
    elif not os.path.isfile(attributes):
        ConcatenateDatasets("genomes/genomes.txt")
    #PanOCT.RunPanOCT(fasta_db, attributes, blast, tags, **kwargs)

    # If enabled, try to fill potential gaps in syntenic clusters within pangenome using BLAST+ data.
    if gaps:
        logging.info("Master: Running gap filling method.")
        PanOCT.FillGaps(blast, "./matchtable.txt", fasta_db, "./panoct_tags.txt")
        PanOCT.PanOCTOutputHandler()
        PanOCT.GenerateClusterFASTAs("genomes/genomes.txt", gaps)
    else:
        pass
        PanOCT.PanOCTOutputHandler()
        PanOCT.GenerateClusterFASTAs("genomes/genomes.txt")


def IPSHandler(cores=None):
    """
    Run InterProScan annotation of pangenome dataset. Note, this only works on Linux and won't run otherwise.
    """
    GO.RunInterProScan("/gm_pred/sets/allprot.db", cores)


def GOHandler(refined=False):
    """
    Run GO-slim enrichment analysis on pangenome datasets using GOATools.
    """
    # Make GO folder, unless one already exists.
    GO.MakeWorkingDirs()

    # Generate dictionary for IPS annotation data.
    annos = GO.GenerateAnnoDict("ips.output.tsv")

    # Generate GO associations and populations files.
    GO.GenerateAssociations(annos)
    if refined:
        GO.GeneratePopulations(annos, "./panoct/refined_matchtable.txt")
    else:
        GO.GeneratePopulations(annos, "./panoct/matchtable.txt")

    # Run map_to_slim.py from GOATools.
    GO.GenerateSlimData("go/associations.txt", "go.obo", "goslim_generic.obo")

    # Run enrichment analysis of core and accessory genomes against full pangenome dataset.
    GO.CoreEnrichment("go.obo", "go/core_pop.txt", "go/full_pop.txt", "go/pangenome_slim.txt")
    GO.AccessoryEnrichment("go.obo", "go/acc_pop.txt", "go/full_pop.txt", "go/pangenome_slim.txt")


def PAMLHandler(ml_path, yn_path, refine=False):
    """
    Run Yn00 on core and accessory gene model clusters.
    """
    if refine:
        clusters = glob("./panoct/clusters/refined/core/fna/Core*.fna") + glob("./panoct/clusters/refined/acc/fna/Acc*.fna")
    else:
        clusters = glob("./panoct/clusters/core/fna/Core*.fna") + glob("./panoct/clusters/acc/fna/Acc*.fna")
    for cluster in clusters:
        trans_seqs = PAML.TranslateCDS(cluster)
        prot_alignment = PAML.MUSCLEAlign(ml_path, trans_seqs)
        nucl_alignment = PAML.PutGaps(prot_alignment, cluster)
        PAML.RunYn00(yn_path, nucl_alignment)

    PAML.SummarizeYn00(refine)


def KaryoploteRHandler(refined=False):
    """
    Generates chromosomal plots of core and accessory gene models for each genome in a dataset, similar to
    the Ruby program PhenoGram but with way less overhead.
    """
    # Make lengths and karyotypes files (don't think these can be passed as objects to R without a lot of effort).
    Karyotype.GenerateContigLengths("./genomes")
    if refined:
        Karyotype.GenerateKaryotypeFiles("./gm_pred/sets/allatt.db", "./panoct/refined_matchtable.txt")
    else:
        Karyotype.GenerateKaryotypeFiles("./gm_pred/sets/allatt.db", "./panoct/matchtable.txt")

    # Pass required files to KaryPloteR and run R script.
    Karyotype.KaryoPloteR("./panoct_tags.txt", "./karyotypes.txt", "./genomes/lengths.txt")


def SizeVizHandler(refined=False):
    """
    Generates bar chart plot of syntenic ortholog cluster sizes in a pangenome dataset, counts observed number
    of clusters (i.e. observed pangenome size, N) and uses that to estimate the predicted number of
    syntenic clusters by the Chao lower bound estimate method (Eng or N-hat). Also generates ring chart
    for pangenome size.
    """
    # If refined pangenome dataset has been made, use that as the basis for the bar chart and Chao estimates.
    # Also generate ring charts at this point too.
    if refined:
        Size.GenerateRingChart("./panoct/refined_matchtable.txt")
        Size.GenerateSizeNumbers("./panoct/refined_matchtable.txt")
    else:
        Size.GenerateRingChart("./panoct/matchtable.txt")
        Size.GenerateSizeNumbers("./panoct/matchtable.txt")

    # Pass required files to BarChart.R and run script.
    Size.GenerateBarChart("./cluster_sizes.txt")


def UpSetRHandler(refined=False):
    """
    Generates UpSetR plot for the distribution of syntenic orthologs within a pangenome dataset, similar to a Venn
    or Euler diagram but capable of >7 input sets.
    """
    if refined:
        UpSet.UpSetR("./panoct_tags.txt", "./panoct/refined_matchtable.txt")
    else:
        UpSet.UpSetR("./panoct_tags.txt", "./panoct/matchtable.txt")


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

    # Add arguments for optional quality control analyses.
    ap.add_argument("--qc", action="store_true", help="Perform quality check on predicted gene model sets.")
    ap.add_argument("--busco", action="store_true", help="Perform BUSCO analysis on predicted gene model sets.")

    # Add argument for skipping all-vs.-all BLASTp step (usually faster to generate data elsewhere).
    ap.add_argument("--no_blast", action="store_true", help="Skip all-vs.-all BLASTp step for PanOCT.")

    # Add argument for skipping PanOCT analysis (mostly for debugging purposes).
    ap.add_argument("--no_panoct", action="store_true", help="Skip PanOCT analysis.")

    # Add argument for gap filling in PanOCT-dervied pangenome.
    ap.add_argument("--refine", action="store_true", help="Attempt to fill potential gaps in syntenic clusters.")

    # Add arguments for annotation and GO-enrichment analysis.
    ap.add_argument("--ips", action="store_true", help="Perform InterProScan analysis of gene model sets. NOTE:"
                                                       " Do not enable this option on non-Linux operating systems,"
                                                       " InterProScan is not supported on these systems.")
    ap.add_argument("--go", action="store_true", help="Perform GO-slim enrichment analysis using GOATools.")

    # Add argument for selection analysis using yn00.
    ap.add_argument("--yn00", action="store_true", help="Perform selection analysis on core and accessory gene "
                                                        "families using yn00.")

    # Add argument to produce all R plots.
    ap.add_argument("--plots", action="store_true", help="Generate all downstream plots (karyotype, cluster size, "
                                                         "ring chart, UpSet, &c).")

    # Add argument to produce karyotype plots.
    ap.add_argument("--karyo", action="store_true", help="Generate karyotype plots for all genomes in a "
                                                         "database based on PanOCT results.")

    # Add argument to produce ring chart of pangenome component size and bar charts with
    # ortholog cluster sizes and Chao (1987) estimates of true pangenome size.
    ap.add_argument("--size", action="store_true", help="Generate ring and bar charts of pangenome complement, "
                                                        "observed and predicted sizes.")

    # Add argument to produce UpSet plot of distribution of syntenic orthologs within accessory genome.
    ap.add_argument("--upset", action="store_true", help="Generate UpSet plot of distribution of syntenic orthologs "
                                                         "within accessory genome of a pangenome dataset.")

    #ap.add_argument("--order", action="store_true", help="Generate circos plot of cluster order within pangenome "
    #                                                     "dataset (implies --karyo).")

    # Add mandatory positional argument for path to config file (default will be the .ini file in /src).
    ap.add_argument("CONFIG_FILE", action="store", nargs="?", help="Path to PanGLOSS configuration file.",
                    default=os.path.dirname(os.path.realpath(sys.argv[0])) + "/config.ini")

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

    # Read config file (either the default or what's supplied at the command-line).
    logging.info("Master: Read config file.")
    if not ap.CONFIG_FILE:
        ap.CONFIG_FILE = os.path.dirname(os.path.realpath(sys.argv[0])) + "/config.ini"
    cp.read(ap.CONFIG_FILE)

    # Get paths for prediction dependencies.
    for arg in cp.items("Gene_prediction_dependencies"):
        if arg[0] == "exonerate_path":
            ex_path = arg[1]
        if arg[0] == "genemark_path":
            gm_path = arg[1]
        if arg[0] == "transdecoder_predict_path":
            tp_path = arg[1]
        if arg[0] == "transdecoder_longorfs_path":
            tl_path = arg[1]

    # Get paths for QC dependencies.
    for arg in cp.items("Quality_control_dependencies"):
        if arg[0] == "busco_path":
            bu_path = arg[1]
        if arg[0] == "busco_lineage_path":
            bl_path = arg[1]

    # Get paths for downstream analysis dependencies.
    for arg in cp.items("Analysis_dependencies"):
        if arg[0] == "muscle_path":
            ml_path = arg[1]
        if arg[0] == "yn00_path":
            yn_path = arg[1]
        if arg[1] == "ips_path":
            ip_path = arg[1]

    # Unless disabled, parse arguments for PanGuess and run gene model prediction.
    if ap.pred or ap.pred_only:
        panguess_args = [ex_path, gm_path, tp_path, tl_path]
        logging.info("Master: Performing gene prediction steps using PanGuess.")
        for arg in cp.items("Gene_model_prediction"):
            panguess_args.append(arg[1])
        if ap.no_exonerate:
            panguess_args.append(True)
        PanGuessHandler(*panguess_args)
        logging.info("Master: Gene prediction finished.")

        # If enabled, check gene sets against user-provided sets of dubious genes, or transposable elements, &c.
        if ap.qc:
            logging.info("Master: Performing gene model QC using QualityCheck.")
            qc_args = [[i for i in glob("./gm_pred/sets/*.faa")]]
            for arg in cp.items("Quality_control"):
                if arg[1]:
                    qc_args.append(arg[1])
            QualityCheckHandler(*qc_args)
            logging.info("Master: QC analysis finished.")
        else:
            logging.info("Master: Skipped gene model QC (--qc not enabled).")

        # If enabled, run BUSCO completedness analysis of all gene model sets.
        if ap.busco:
            logging.info("Master: Performing BUSCO analysis of gene model sets.")
            busco_args = [bu_path, bl_path]
            busco_args = busco_args + [[i for i in glob("./gm_pred/sets/*.faa")]]
            BUSCOHandler(*busco_args)

        # Allow program to finish after gene prediction and (optionally) QC/BUSCO if --pred_only is enabled.
        if ap.pred_only:
            logging.info("Master: Finishing PanGLOSS (--pred_only enabled). To run remaining steps with your own "
                         "all-vs.-all BLAST data, run PanGLOSS with the --no_pred and --no_blast flags.")
            logging.info("Master: PanGLOSS finished in {0}.".format(str(datetime.now() - start_time)))
            sys.exit(0)
    else:
        logging.info("Master: Skipped gene prediction steps (--nopred enabled).")

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

    # Run PanOCT on full dataset, unless --no_panoct is enabled.
    if not ap.no_panoct:
        logging.info("Master: Performing PanOCT analysis of dataset.")
        panoct_default_args = []
        panoct_extra_args = []
        for arg in cp.items("PanOCT_settings"):
            if arg[1]:
                panoct_default_args.append(arg[1])
        if ap.refine:
            panoct_default_args.append(True)
        if panoct_extra_args:
            pass
        else:
            PanOCTHandler(*panoct_default_args)
    else:
        logging.info("Master: Skipping PanOCT analysis (--no_panoct enabled).")

    # If enabled, run InterProScan analysis on entire dataset.
    if ap.ips:
        if sys.platform != "linux":
            print "InterProScan is not supported on non-Linux operating systems. Cannot run InterProScan analysis."
            print "See https://github.com/ebi-pf-team/interproscan/wiki for more information."
            pass
        else:
            IPSHandler()

    # If enabled, run GO-slim enrichment analysis on core and accessory datasets using GOATools.
    if ap.go:
        GOHandler(ap.refine)
        pass

    # If enabled, run selection analysis using yn00.
    if ap.yn00:
        logging.info("Master: Performing selection analysis using yn00.")
        PAMLHandler(ml_path, yn_path, ap.refine)

    # If enabled, enable all plot arguments.
    if ap.plots:
        ap.karyo = True
        ap.size = True
        ap.upset = True

    # If enabled, generate karyotype plots for all strain genomes in pangenome dataset.
    if ap.karyo:
        logging.info("Master: Generating karyotype plots for all genomes in dataset.")
        KaryoploteRHandler(ap.refine)

    # If enabled, generate bar charts and Chao estimate of pangenome size.
    if ap.size:
        logging.info("Master: Generating size plots.")
        SizeVizHandler(ap.refine)

    # If enabled, generate UpSet plot of accessory genome.
    if ap.upset:
        logging.info("Master: Generating UpSet accessory genome distribution plot.")
        UpSetRHandler(ap.refine)


if __name__ == "__main__":
    main()
