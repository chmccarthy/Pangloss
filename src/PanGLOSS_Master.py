from argparse import ArgumentParser
from Bio import SearchIO, SeqIO
from ConfigParser import SafeConfigParser
from PanGLOSS import PanGuess

def PanGuessHandler(genomelist, workdir, ref):
    """
    Runs PanGuess from master script.
    """
    PanGuess.MakeWorkingDir(workdir)
    PanGuess.BuildRefSet(workdir, ref)
    

### Parser functions. ###

def CmdlineParser():
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
    ap = CmdlineParser()
    cp = ConfigFileParser()
    
    # Create argument lists.
    panguess_args = []
    pangloss_args = []
    
    
    
    # Read config file.
    cp.read(ap.config)
    
    # Generate arguments for PanGuess (genome list, working directory, reference set).
    for arg in cp.items("Gene_model_prediction"):
        panguess_args.append(arg[1])
    
    # Run PanGuess
    PanGuessHandler(*panguess_args)

if __name__ == "__main__":
    main()