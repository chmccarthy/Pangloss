from argparse import ArgumentParser
from ConfigParser import SafeConfigParser
from glob import glob
from PanGLOSS import PanGuess, QualityCheck


def PanGuessHandler(genomelist, workdir, ref, exon_cov, gm_branch, td_potenial, td_len):
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

    Find some way to incorporate exon_cov and td_potential into final product!
    """
    # Generate list of genomes from user-provided genome list file.
    print "Parsing list of genomes...\t",
    genomes = [line.strip("\n") for line in open(genomelist)]
    print "OK."
    
    # Create working directory and split reference protein set if not already present.
    print "Creating working directory and building reference protein set...\t",
    PanGuess.MakeWorkingDir(workdir)
    PanGuess.BuildRefSet(workdir, ref)
    print "OK."
    
    # Loop over each genome and carry out gene model prediction.
    for genome in genomes:
        # Make tag from genome name.
        print "Running gene model prediction for {0}.".format(genome)
        tag = genome.split(".")[0].split("/")[1]
        
        # Run prediction using Exonerate.
        print "Performing gene model prediction for {0} using Exonerate...\t".format(tag),
        cmds = PanGuess.BuildExonerateCmds(workdir, genome)
        exonerate_genes = PanGuess.RunExonerate(cmds)
        print "OK."
        
        # Order gene models predicted via Exonerate by Contig ID: Location.
        print "Sorting predicted Exonerate gene models by genomic location...\t",
        exonerate_genes.sort(key=lambda x: (x.contig_id, x.locs[0]))
        print "OK."
        
        # Extract genomic attributes from Exonerate gene model set.
        print "Extracting genomic attributes from Exonerate gene models...\t",
        exonerate_attributes = PanGuess.GetExonerateAttributes(exonerate_genes, tag)
        print "OK."
        
        # Run prediction using GeneMark-ES.
        print "Running gene model prediction for {0} using GeneMark-ES...\t".format(genome),
        genemark_gtf = PanGuess.RunGeneMark(genome, gm_branch)
        print "OK."
        
        # Convert GeneMark-ES GTF file into a more PanOCT-compatible version.
        print "Converting GeneMark-ES GTF file to attributes...\t",
        genemark_attributes = PanGuess.GeneMarkGTFConverter(genemark_gtf, tag)
        print "OK."
        
        # Merge unique gene model calls between Exonerate and GeneMark-ES.
        print "Merging unique gene calls...\t",
        merged_attributes = PanGuess.MergeAttributes(exonerate_attributes, genemark_attributes)
        print "OK."
        
        # Clean up GeneMark-ES files and folders.
        print "Tidying up GeneMark-ES temporary files and folders...\t",
        PanGuess.MoveGeneMarkFiles(workdir, genome)
        print "OK."
        
        # Extract NCRs into list.
        print "Extracting non-coding regions of genome...\t",
        noncoding = PanGuess.ExtractNCR(merged_attributes, genome)
        print "OK."
        
        # Run TransDecoder on NCRs.
        print "Running TransDecoder on non-coding regions...\t",
        tdir = PanGuess.RunTransDecoder(noncoding, workdir, genome, td_len)
        print "OK."
        
        # Move TransDecoder files.
        print "Tidying up TransDecoder temporary files and folders...\t",
        PanGuess.MoveTransDecoderFiles(tdir)
        print "OK."
        
        # Extract TransDecoder attributes.
        print "Converting TransDecoder GTF file to attributes...\t",
        trans_attributes = PanGuess.TransDecoderGTFToAttributes(tdir, tag)
        print "OK."
        
        # Merge TransDecoder calls into the Exonerate + GeneMark-ES set.
        print "Merging remaining gene calls...\t",
        full_attributes = PanGuess.MergeAttributes(merged_attributes, trans_attributes)
        print "OK."
        
        # Write out gene set, protein set and attributes set.
        print "Writing gene calls and gene attributes...\t",
        PanGuess.ConstructGeneModelSets(full_attributes, exonerate_genes, workdir, genome, tag)
        print "OK."
        
        # Compress temporary folders and finish up.
        print "Compressing temporary GeneMark-ES and TransDecoder folders...\t",
        PanGuess.TarballGenePredictionDirs(workdir, genome)
        print "OK."
        print "Finished gene model prediction for {0}.".format(genome)


def QualityCheckHandler(gene_sets, queries):
    """
    Search a user-provided set of genes of dubious-quality (i.e. pseudogenes, transposable elements or
    transposons &c.) against predicted gene model sets and filter out sufficiently similar genes in the latter.

    Arguments:
        gene_sets   = List of strains in analysis (easy access to all files associated with a strain).
        queries     = Set of genes (protein sequences, in fact) to search against all gene model sets.
    """
    QualityCheck.BuildMakeBLASTDBs(gene_sets)
    results = QualityCheck.QCBLAST(queries, gene_sets)
    for result in results:
        for q in result:
            print q
    pass


def BUSCOHandler():
    pass


### Parser functions. ###

def CmdLineParser():
    """
    Create and return a configuration file parser.
    """
    ap = ArgumentParser(description="Pan-genome analysis of microbial eukaryotes.")
    ap.add_argument("config", help="Path to PanGLOSS configuration file.")
    ap.add_argument("--no_pred", action="store_true", help="Skip gene model prediction pipeline.")
    ap.add_argument("--qc", action="store_true", help="Perform quality check on predicted gene model sets.")
    ap.add_argument("--busco", action="store_true", help="Perform BUSCO analysis on predicted gene model sets.")
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
    # Create parser objects.
    ap = CmdLineParser()
    cp = ConfigFileParser()
    
    # Create argument lists.
    panguess_args = []
    pangloss_args = []

    # Read config file.
    cp.read("config.txt")

    # Unless disabled, parse arguments for PanGuess and run gene model prediction.
    if not ap.no_pred:
        for arg in cp.items("Gene_model_prediction"):
            panguess_args.append(arg[1])
        PanGuessHandler(*panguess_args)

    # If enabled, check gene sets against user-provided sets of dubious genes, or transposable elements.
    if ap.qc:
        gene_sets = [i.split(".")[0] for i in glob("*.faa")]
        for arg in cp.items("Quality_check"):
            queries = arg[1]
        QualityCheckHandler(gene_sets, queries)


if __name__ == "__main__":
    main()

