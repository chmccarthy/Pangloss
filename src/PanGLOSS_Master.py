from argparse import ArgumentParser
from Bio import SearchIO, SeqIO
from ConfigParser import SafeConfigParser
from PanGLOSS import PanGuess
from csv import reader

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
        
        # Extract genomic attributes from Exonerate gene model set (easier to do at
        # this point rather than later).
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
        
        # Merge unique gene model calls between the two different methods.
        print "Merging unique gene calls...\t",
        merged_attributes = PanGuess.MergeAttributes(tag, exonerate_attributes, genemark_attributes)
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
        tdir = PanGuess.RunTransDecoder(noncoding, workdir, genome)
        print "OK."
        
        # Move TransDecoder files.
        print "Tidying up TransDecoder temporary files and folders...\t",
        PanGuess.MoveTransDecoderFiles(tdir)
        print "OK."
        
        # Extract TransDecoder attributes.
        print "Converting TransDecoder GTF file to attributes...\t"
        trans_attributes = PanGuess.TransDecoderGTFToAttributes(tdir, tag)
        print "OK."
        
        # Merge unique gene model calls between the two different methods.
        print "Merging remaining gene calls...\t",
        full_attributes = PanGuess.MergeAttributes(tag, merged_attributes, trans_attributes)
        print "OK."
        
        print "Writing gene calls and gene attributes...\t"
        PanGuess.ConstructGeneModelSets(full_attributes, exonerate_genes, workdir, genome, tag)
        print "OK."


def QualityCheck():
    pass


### Parser functions. ###

def CmdLineParser():
    """
    Create and return a configuration file parser.
    """
    ap = ArgumentParser(description="Pan-genome analysis of microbial eukaryotes.")
    ap.add_argument("config", help="Path to PanGLOSS configuration file.")
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
    
    # Generate arguments for PanGuess (genome list, working directory, reference set).
    for arg in cp.items("Gene_model_prediction"):
        panguess_args.append(arg[1])
    
    # Run PanGuess, unless disabled.
    PanGuessHandler(*panguess_args)

if __name__ == "__main__":
    main()