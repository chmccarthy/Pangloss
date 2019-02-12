#!/usr/bin/env perl
#Copyright (C) 2011-2015 The J. Craig Venter Institute (JCVI).  All rights reserved
#Written by Derrick E. Fouts, Ph.D. and Granger Sutton, Ph.D.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#Revision notes
# Changed shebang line from /usr/local/bin/perl to "/usr/bin/env perl" for increased portablity
# Using cwd() instead of $ENV{'PWD'} for increased portability.
# Made warning for case when best BLAST score > than self score show only in debug mode
# Increased input and output file documentation in README.txt file
# Added options T, C, W, S, U. B, c and e
# Fixed bugs (refer to the svn code repository)
# Match table and other revisions 01/07/2004
# Added NCBI -m9 input option 12/08/2010
# Moved genome identification from feat/locus_name-tag to the gene_att file 12/08/2010
my $commandline = join (" ", @ARGV);
my $prog = $0;
$prog =~ s/.*\///;

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);
getopts ('TdhDb:p:t:f:i:I:E:L:M:G:H:V:P:Q:N:g:A:C:F:a:s:W:S:U:B:c:e:R:');#GGS added F, a, s, options F option for frameshift and a option for amount of missing amino acids to still be full length and s option for an evidence threshold for frameshifts
our ($opt_T,$opt_d,$opt_h,$opt_D,$opt_b,$opt_p,$opt_t,$opt_f,$opt_i,$opt_I,$opt_E,$opt_L,$opt_M,$opt_G,$opt_H,$opt_V,$opt_P,$opt_Q,$opt_N,$opt_g,$opt_A,$opt_C,$opt_F,$opt_a,$opt_s,$opt_W,$opt_S,$opt_U, $opt_B, $opt_c, $opt_e, $opt_R);

my ($print_cluster_numbers,$compute_centroids,$basedir,$btabpath,$pep_path,$att_file,$pep_file,$btabfile,$tagfile,$percentid,$fspercentid,$evalue,$min_hit_length,$microarray,$histogramfile,$hitsfile,$vennfile,$normalizefile,$writeparalogs,$writeclusterweights,$DEBUG,$frameshiftfile,$frameshift_length_ratio,$max_missing_aa,$frames_link_thresh,$ANCHOR_window_size,$ANCHOR_cutoff,$CGN_window_size,$strict_orthos,$fragmentfile, $pairwisefile, $compute_adjacency, $cluster_input_file, $process_clusters);#GGS added frameshift variables
#GGS $frameshiftfile is a boolean for detecting and outputting frameshifted or truncated protein fragments true if -F
#GGS -F $frameshift_length_ratio is a value between 1 and 2 (default 1.33) used to determine that we are detecting adjacent protein fragments rather than adjacent paralogs
#GGS -a $max_missing_aa is how many amino acids at either end of a blast match can be missing and still be considered full length
#GGS -s $frames_link_thresh is a threshold for number of blast matches indicating a frameshift to believe it
#GGS -G yes outputs a histograms file
#set defaults
$microarray = 0;
$hitsfile = 0;
$histogramfile = 0;
$strict_orthos = 1;
$normalizefile = 0;
$writeparalogs = 1;
$fragmentfile = 1;
$pairwisefile = 1;
$writeclusterweights = 1;
$vennfile = 1;
$compute_adjacency = 0;
$compute_centroids = 1;
$ANCHOR_window_size = 2;
$ANCHOR_cutoff = (2 * $ANCHOR_window_size) + 1;
$cluster_input_file = "";
$process_clusters = 0;

my @cluster_levels = ();
my @BSR_cluster_levels = ();

## use boolean logic:  TRUE = 1, FALSE = 0

my $version = "3.23";
if ($opt_d) {
    $compute_centroids = 1;
} else {
    $compute_centroids = 1; #default - always compute centroids
}
if ($opt_T) {
    $print_cluster_numbers = 1;
} else {
    $print_cluster_numbers = 0; #default
}
if ($opt_D) {$DEBUG = 1;} else { $DEBUG = 0; } # Debug mode is off as default.
#GGS process strict ortholog option
if ($opt_S) {
    if (($opt_S eq "Y") or ($opt_S eq "y") or ($opt_S eq "Yes") or ($opt_S eq "yes") or ($opt_S eq "YES")) {$strict_orthos = 1;} 
    elsif (($opt_S eq "N") or ($opt_S eq "n") or ($opt_S eq "No") or ($opt_S eq "no") or ($opt_S eq "NO")) {$strict_orthos = 0;}
    elsif (($opt_S eq "M") or ($opt_S eq "m") or ($opt_S eq "More") or ($opt_S eq "more") or ($opt_S eq "MORE")) {$strict_orthos = 2;}
    else { $strict_orthos = 1; }  # Want to use strict criteria for orthologs  0 = NO, 1 = YES [DEFAULT = YES]
}
#GGS process CGN window size option
if ($opt_W) {
    if (($opt_W >= 1) && ($opt_W <= 20)) {#GGS ratioCGN window size must be strictly between 1 and 20
	$CGN_window_size = int $opt_W;
    } else {
	&option_help;
    }
} else {
    $CGN_window_size = 5; #default CGN window size
}
#GGS process frameshift option
if ($opt_F) {
    $frameshiftfile = 1;
    if (($opt_F < 2) && ($opt_F > 1)) {#GGS ratio must be strictly between 1 and 2
	$frameshift_length_ratio = $opt_F;
    } else {
	&option_help;
    }
} else {
    $frameshiftfile = 0;
}
#GGS set max_missing_aa
if ($opt_a) {
    if (($opt_a > 0) && ($opt_a <= 100)) {#GGS must be between 0 and 100
	$max_missing_aa = $opt_a;
    } else {
	&option_help;
    }
} else {
    $max_missing_aa = 20;#default
}
#GGS set frameshift evidence threshold
if ($opt_s) {
    if ($opt_s > 0) {#GGS must be > 0
	$frames_link_thresh = $opt_s;
    } else {
	&option_help;
    }
} else {
    $frames_link_thresh = 1;#default
}
if ($opt_h) { &option_help; } # quit with help menu
if ($opt_b) {$basedir = $opt_b;} else { $basedir = $ENV{'PWD'}; } # if no value for option b (base or working directory) set it to current directory
if ($opt_p) {$btabpath = $opt_p;} else { $btabpath = $basedir; } # if no value given for path to btab file, make default = base directory
if ($opt_Q) {$pep_path = $opt_Q;} else { $pep_path = $basedir; } # if no value given for path to peptide file, make default = base directory
if (($opt_g) && (-s "$opt_g")) {$att_file = $opt_g;} else { print STDERR "Error with -g\n"; &option_help; } # if no value for option g (name of gene_att file), quit with help menu
if (($opt_P) && (-s "$pep_path/$opt_P")) {$pep_file = $opt_P;} else { print STDERR "Error with -P $pep_path/$opt_P\n"; &option_help; } # if no value for option P (pepfile), quit with help menu
if (($opt_f) && (-s "$opt_f")) {$tagfile = $opt_f;} else { print STDERR "Error with -f $opt_f\n"; &option_help; }  # must supply the file containing the list of unique genome names to analyze
if ($opt_R) {
    if ($opt_t) {
	print STDERR "Blast tabular output file: $opt_t is being ignored since previously computed clusters are being used (-R option)\n";
    }
    if (-s "$opt_R") {
	$cluster_input_file = $opt_R;
	$process_clusters = 1;
	$btabfile = ""; #so outputting the parameters file doesn't complain
    } else { #if supplied the -R argument must be a file of PanOCT formatted clusters
	print STDERR "Error with -R $opt_R\n";
	&option_help;
    }
} else {
    $process_clusters = 0;
    #only need blast data if we are not just processing clusters
    if (($opt_t) && (-s "$btabpath/$opt_t")) {$btabfile = $opt_t;} else { print STDERR "Error with -t $btabpath/$opt_t\n"; &option_help; } # if no value for option t (name of btab file), quit with help menu
}
if ($opt_i) {$percentid = $opt_i;} else { $percentid = 35.0; } # the minimum cutoff to use blast scores for possible matches
if ($opt_I) {$fspercentid = $opt_I;} else { $fspercentid = 35.0; } # the minimum cutoff to use blast scores for frameshift detection
if ($opt_E) {$evalue = $opt_E;} else { $evalue = 0.00001; } # if no E-value cut-off given, make default 0.00001
if ($opt_L) {$min_hit_length = $opt_L;} else { $min_hit_length = 1; }
if ($opt_M) {
    if (($opt_M eq "Y") or ($opt_M eq "y") or ($opt_M eq "Yes") or ($opt_M eq "yes") or ($opt_M eq "YES"))  {$microarray = 1;} 
    elsif (($opt_M eq "N") or ($opt_M eq "n") or ($opt_M eq "No") or ($opt_M eq "no") or ($opt_M eq "NO")) {$microarray = 0;}
    else { $microarray = 0; }  # Want to create microarray-like data of normalized BLAST scores  0 = NO, 1 = YES [DEFAULT = NO]
}
if ($opt_G) {
    if (($opt_G eq "Y") or ($opt_G eq "y") or ($opt_G eq "Yes") or ($opt_G eq "yes") or ($opt_G eq "YES")) {$histogramfile = 1;} 
    elsif (($opt_G eq "N") or ($opt_G eq "n") or ($opt_G eq "No") or ($opt_G eq "no") or ($opt_G eq "NO")) {$histogramfile = 0;}
    else { $histogramfile = 0; }  # Want to create histograms  0 = NO, 1 = YES [DEFAULT = NO]
}
if ($opt_H) {
    if (($opt_H eq "Y") or ($opt_H eq "y") or ($opt_H eq "Yes") or ($opt_H eq "yes") or ($opt_H eq "YES")) {$hitsfile = 1;} 
    elsif (($opt_H eq "N") or ($opt_H eq "n") or ($opt_H eq "No") or ($opt_H eq "no") or ($opt_H eq "NO")) {$hitsfile = 0;}
    else { $hitsfile = 0; }  # Want to create table of hits  0 = NO, 1 = YES [DEFAULT = NO]
}
if ($opt_V) {
    if (($opt_V eq "Y") or ($opt_V eq "y") or ($opt_V eq "Yes") or ($opt_V eq "yes") or ($opt_V eq "YES")) {$vennfile = 1;}
    elsif (($opt_V eq "N") or ($opt_V eq "n") or ($opt_V eq "No") or ($opt_V eq "no") or ($opt_V eq "NO")) {$vennfile = 0;}
    else { $vennfile = 1; } # Want to create match table file?  0 = NO, 1 = YES [DEFAULT = YES]
}
if ($opt_N) {
    if (($opt_N eq "Y") or ($opt_N eq "y") or ($opt_N eq "Yes") or ($opt_N eq "yes") or ($opt_N eq "YES")) {$normalizefile = 1;}
    elsif (($opt_N eq "N") or ($opt_N eq "n") or ($opt_N eq "No") or ($opt_N eq "no") or ($opt_N eq "NO")) {$normalizefile = 0;}
    else { $normalizefile = 0; } # Want to create normalized BLAST score file?  0 = NO, 1 = YES [DEFAULT = NO]
}
if ($opt_U) {
    if (($opt_U eq "Y") or ($opt_U eq "y") or ($opt_U eq "Yes") or ($opt_U eq "yes") or ($opt_U eq "YES")) {$fragmentfile = 1;}
    elsif (($opt_U eq "N") or ($opt_U eq "n") or ($opt_U eq "No") or ($opt_U eq "no") or ($opt_U eq "NO")) {$fragmentfile = 0;}
    else { $fragmentfile = 1; } # Want to create list of possible protein fragments and fusions?  0 = NO, 1 = YES [DEFAULT = YES]
}
if ($opt_B) {
    if (($opt_B eq "Y") or ($opt_B eq "y") or ($opt_B eq "Yes") or ($opt_B eq "yes") or ($opt_B eq "YES")) {$pairwisefile = 1;}
    elsif (($opt_B eq "N") or ($opt_B eq "n") or ($opt_B eq "No") or ($opt_B eq "no") or ($opt_B eq "NO")) {$pairwisefile = 0;}
    else { $pairwisefile = 1; } # Want to create list of pairwise similarity scores?  0 = NO, 1 = YES [DEFAULT = YES]
}
if  ( $opt_A ) { 
    if(($opt_A eq "Y") or ($opt_A eq "y") or ($opt_A eq "Yes") or ($opt_A eq "yes") or ($opt_A eq "YES")) {$writeparalogs = 1;}
    elsif (($opt_A eq "N") or ($opt_A eq "n") or ($opt_A eq "No") or ($opt_A eq "no") or ($opt_A eq "NO")) {$writeparalogs = 0;}
    else { $writeparalogs = 1; } # Want to create a list of paralogs?  0 = NO, 1 = YES [DEFAULT = YES]
}
if  ( $opt_C ) { 
    if(($opt_C eq "Y") or ($opt_C eq "y") or ($opt_C eq "Yes") or ($opt_C eq "yes") or ($opt_C eq "YES")) {$writeclusterweights = 1;}
    elsif (($opt_C eq "N") or ($opt_C eq "n") or ($opt_C eq "No") or ($opt_C eq "no") or ($opt_C eq "NO")) {$writeclusterweights = 0;}
    else { $writeclusterweights = 1; } # Want to create a list of cluster weights?  0 = NO, 1 = YES [DEFAULT = YES]
}
#process percent of representation in cluster for adjacency matrix into cluster_levels array
if ($opt_c) {
    $compute_adjacency = 1;
    @cluster_levels = split(',', $opt_c);
    foreach my $level (@cluster_levels) {
	if ($DEBUG) {
	    print STDERR "calc_adjacency level $level\n";
	}
	if (!(looks_like_number($level))) {
	    die ("ERROR: $level is not a number in list of adjacency matrix cluster representation values (-c)!\n");
	}
	if (($level > 100) || ($level < 0)) {
	    die ("ERROR: $level is not between 0-100 in list of adjacency matrix cluster representation values (-c)!\n");
	}
    }
}
#process percent of representation in cluster for BSR distance matrices into BSR_cluster_levels array
if ($opt_e) {
    @BSR_cluster_levels = split(',', $opt_c);
    foreach my $level (@BSR_cluster_levels) {
	if ($DEBUG) {
	    print STDERR "BSR level $level\n";
	}
	if (!(looks_like_number($level))) {
	    die ("ERROR: $level is not a number in list of BSR matrix cluster representation values (-e)!\n");
	}
	if (($level > 100) || ($level < 0)) {
	    die ("ERROR: $level is not between 0-100 in list of BSR matrix cluster representation values (-e)!\n");
	}
    }
} else {
    $BSR_cluster_levels[0] = 0;
    $BSR_cluster_levels[1] = 90;
    $BSR_cluster_levels[2] = 100;
}

my %query_evalue_cutoff = ();# Key = query id, Value = query specific evalue cutoff for low scoring queries
my @paralog_cutoff = ();    # Index1 = genome tag index, Index2 = genome tag index, Value = [0,1] normalized blast score cutoff for paralogs between these two genomes
my @ave_per_id = ();        # Index1 = genome tag index, Index2 = genome tag index, Value = [0,100] average percent identity for high quality orthologs between these two genomes
my @mean_max_BSR = ();      # Index1 = genome tag index, Index2 = genome tag index, Value = [0,1] mean max BSR for high quality orthologs between these two genomes
my %genome_hash = ();       # Key1 = genome tag, Key2 = asmbl_id, Value = array of protein ids sorted by end5 and then end3
my %genome_hash_context = ();# Key1 = genome tag, Key2 = asmbl_id, Key3 = positions in genome_hash that are only for context, Value = 1 just to define the key3
my %feat_hash = ();         # Key = feat_name Key2 = struct members with their values
my %relationship_hash = (); # Key1 = Query protein, Key2 = Subject protein, Key3 = struct members with their values
my %Qbytaghash = ();        # Key1 = Query protein, Key2 = subject genome tag, Key3 = Subject protein, Key4 = same as relationship_hash Key3 (actually is the same pointer)
my %BestTagMatch = ();      # Key1 = Query protein, Key2 = subject genome tag, Value = bit score of best match for this genome
my %SecondBestTagMatch = ();# Key1 = Query protein, Key2 = subject genome tag, Value = bit score of best match for this genome
my %orf_counter = ();       # Key = genome tag, Value = orf counts
my %Tagbyfeatnamehash = (); # Key1 = genome tag, Key2 = feat_name (query)
my %TagByPointer = ();      # Key = feat_name-tag, value = array position location
my @tag_array = ();         # array of genome tags
my %FeatnameLookupTag_hash = ();# Key = feat_name, Value = tag.  Needed for input lacking tagged feat_names
my %TagIndex = ();    # Key = genome tag, Value = index in tag_array of associated genome tag
my %AssemblyLookup_hash = ();# Generated to store the asmbl_id (value) of each feat_name (key) so we can lookup the correct array during synteny searches
my %clusters = ();          # Key = feat_name, Value = array of (hashes of feat_name and genome_tag) that are in a cluster with the feat_name key
my %cluster_number = ();    # Key = feat_name, Value = cluster number that feat_name is in
my @cluster_size = ();      # Index = cluster number, Value = size of the cluster
my @cluster_for_number =(); # Index = cluster number, Value = cluster structure (array)
my @centroids = ();          # Index = cluster number, Value = feature id for the centroid of the cluster
my $minSynMatches = 2;      # Hard-code the minimum number of matches to be consider syntenous
my $min_synteny_threshold = 0;# minimum synteny score to trust
my $genome_number;          # number of genomes defined in tagfile
#my $outprefix = "$basedir/panoct_"; #previous naming convention
my $outprefix = "$basedir/";
my $btab_style = ""; # 1 = WU, 0 = NCBI

#%relationship_hash stores information for the blast matches
###key1 = query
###key2 = subject
###value -> id = percent identity
###         eval = e-value
###         score = BLAST bits score
###         BSR = BLAST Score Ratio # added 12/13/2010 so that the true top matches were chosen, not the Blast score
###         best = 1 if best blast match in the subject genome, otherwise 0
###         bibest = 1 if best bidirectional blast match in the subject genome, otherwise 0
###         min_query = smaller match coordinate on the query protein
###         max_query = larger match coordinate on the query protein
###         min_sub = smaller match coordinate on the subject protein
###         max_sub = larger match coordinate on the subject protein
###         clique_top = number in the clique at the top of blast matches
###         clique_all = number in the clique for all blast matches


sub print_parameters { # print the parameters input via command line or set as default

    my ($outfile) = @_;

    print $outfile "$prog $commandline\n";
    print $outfile "version = $version\n";
    print $outfile "microarray = $microarray\n";
    print $outfile "hitsfile = $hitsfile\n";
    print $outfile "histogramfile = $histogramfile\n";
    print $outfile "strict_orthos = $strict_orthos\n";
    print $outfile "normalizefile = $normalizefile\n";
    print $outfile "writeparalogs = $writeparalogs\n";
    print $outfile "fragmentfile = $fragmentfile\n";
    print $outfile "pairwisefile = $pairwisefile\n";
    print $outfile "print_cluster_numbers = $print_cluster_numbers\n";
    print $outfile "vennfile = $vennfile\n";
    print $outfile "compute_centroids = $compute_centroids\n";
    print $outfile "ANCHOR_window_size = $ANCHOR_window_size\n";
    print $outfile "ANCHOR_cutoff = $ANCHOR_cutoff\n";
    print $outfile "DEBUG = $DEBUG\n";
    print $outfile "CGN_window_size = $CGN_window_size\n";
    print $outfile "frameshiftfile = $frameshiftfile\n";
    print $outfile "max_missing_aa = $max_missing_aa\n";
    print $outfile "frames_link_thresh = $frames_link_thresh\n";
    print $outfile "basedir = $basedir\n";
    print $outfile "btabpath = $btabpath\n";
    print $outfile "pep_path = $pep_path\n";
    print $outfile "att_file = $att_file\n";
    print $outfile "pep_file = $pep_file\n";
    print $outfile "btabfile = $btabfile\n";;
    print $outfile "tagfile = $tagfile\n";
    print $outfile "percentid = $percentid\n";
    print $outfile "fspercentid = $fspercentid\n";
    print $outfile "evalue = $evalue\n";
    print $outfile "min_hit_length = $min_hit_length\n";
    print $outfile "writeclusterweights $writeclusterweights\n";
    print $outfile "compute_adjacency $compute_adjacency\n";
    print $outfile "process_clusters $process_clusters\n";
    print $outfile "cluster_input_file $cluster_input_file\n";
    my $cl_levs = join (" ", @cluster_levels);
    print $outfile "compute_adjacency cluster_levels $cl_levs\n";
    $cl_levs = join (" ", @BSR_cluster_levels);
    print $outfile "BSR_cluster_levels $cl_levs\n";
}

sub get_tags {  # obtain list of genomes to compare
   
    my %temp_hash = ();
    $genome_number = 0;     # total number of genomes to be processed

    open (my $infile, "<", "$basedir/$tagfile") || die ("ERROR: can't open file $basedir/$tagfile\n");
    while (<$infile>)  {
	chomp;
	if (length($_) > 0) {
	    my($name,$location,$asmbl_id) = split(/\t/,$_);
	    $name =~ s/\s+$//;

            if (defined $temp_hash{$name})  {
               die ("ERROR:  You have more than one occurance of $name in $basedir/$tagfile!\n");
            } else  {
		$temp_hash{$name} = 1;
		push (@tag_array, $name); # populate the tag_array in the order of the tagfile (1st tag is the reference tag)
		$genome_number++;
		if ($genome_number == 1) {
		    print STDERR "$name [Reference]\n";
		} else {
		    print STDERR "  $name\n";
		}
	    }
	}  
    }
    close($infile);

    my $index = 0;
    foreach my $tag (@tag_array) {
	$TagIndex{$tag} = $index++
    }
}

sub get_protein_info { #formerly get_fasta_headers.  Added calculation of protein length for ncbi data 12/08/10

  my @line = ();
  my $id;
  my $title = "";
  my $sequence = "";
  my $length = "";

  unless (open (PEPFILE, "<$pep_path/$pep_file") )  {
    die ("can't open file $pep_path/$pep_file.\n");
  }
  my ($save_input_separator) = $/;
  $/="\n>";
  while (<PEPFILE>) {
    ($title,$sequence) = /^>?\s*(.*)\n([^>]+)>?/; # split the header line and sequence (very cool)
    @line = split(/\s+/, $title);  # split the scalar $line on space or tab (to separate the identifier from the header and store in array @fasta
    $id = $line[0]; # unique orf identifier is in column 0, com_name is in rest
    $id =~ s/>//;
    $sequence =~ s/[^a-zA-Z]//g; # remove any non-alphabet letters
    $length = length($sequence);
    $feat_hash{$id}->{'header'} = join(' ', @line[1..$#line]); # put the identifier into the hash as the "key" and the header as value "header" (joining all columns after first space/tab)
    $feat_hash{$id}->{'length'} = $length;
    $feat_hash{$id}->{'sequence'} = $sequence;
    #print STDERR "$id ($feat_hash{$id}->{'header'}) = $feat_hash{$id}->{'length'}\n";
    $title = ""; # clear the title for the next round.
    $sequence = ""; #clear out the sequence for the next round.
  }
  $/ = $save_input_separator; # restore the input separator
  close (PEPFILE);
  return;
}


sub get_gene_att {

    my $pos = "";
    my $tag = "";
    my $end5 = "";
    my $end3 = "";
    my $asmbl_id = "";
    my $feat_name = "";
    my $anno = "";
    my $failed = 0;

    unless (open (ATTFILE, "<$basedir/$att_file") )  {
	die ("ERROR: can not open file $basedir/$att_file.\n");
    }
    while (<ATTFILE>) {
	my @att_line = ();
	chomp;
	@att_line = split(/\t/, $_);  # split the scalar $line on tab
	$asmbl_id = $att_line[0];
	if ($asmbl_id eq "") {
	    print STDERR "ERROR: assembly id/contig id must not be empty/null in the gene attribute file\n$_\n";
	    $failed = 1;
	}
	$feat_name = $att_line[1];
	if (defined $feat_hash{$feat_name}->{'full'}) {
	    print STDERR "ERROR: $feat_name appears more than once in the gene attribute file!\n";
	    $failed = 1;
	}
	$end5 = $att_line[2];
	$end3 = $att_line[3];
	$anno = $att_line[4];
	$tag = $att_line[5];
	if (!defined $TagIndex{$tag}) {
	    print STDERR "ERROR: $tag is a genome tag in the gene attribute file but not in the genome tag file!\n";
	    $failed = 1;
	}
	$feat_hash{$feat_name}->{'full'} = 0; # initialize number of full length blast matches to 0
	$feat_hash{$feat_name}->{'fragment'} = 0; #initialize fragment blast matches to 0
	$feat_hash{$feat_name}->{'fusion'} = 0; #initialize fusion blast matches to 0
	$feat_hash{$feat_name}->{'retained'} = 0; #initialize retained frameshift flag to 0/false
	$FeatnameLookupTag_hash{$feat_name} = $tag;
	push (@{ $genome_hash{$tag}{$asmbl_id} }, $feat_name);
	$AssemblyLookup_hash{$feat_name} = $asmbl_id; # new 01/12/2011
	$feat_hash{$feat_name}->{'end5'} = $end5;
	$feat_hash{$feat_name}->{'end3'} = $end3;
	if ($end5 < $end3) {
	    $feat_hash{$feat_name}->{'orient'} = 1;
	} else {
	    $feat_hash{$feat_name}->{'orient'} = -1;
	}
	if (($feat_name !~ /^CONTEXT[0-9]+:/) && !defined $feat_hash{$feat_name}->{'header'}) {
	    print STDERR "ERROR: $feat_name in the gene attribute file ($att_file) was not found in the protein file ($pep_file)!\n";
	    $failed = 1;
	}
	if ($anno) { # check if there is any annotation from gene_att file
	    $feat_hash{$feat_name}->{'header'} = $anno;
	}
	print STDERR "$feat_name $feat_hash{$feat_name}->{'header'} $feat_hash{$feat_name}->{'end5'} $feat_hash{$feat_name}->{'end3'} $feat_hash{$feat_name}->{'length'} $feat_hash{$feat_name}->{'orient'} $tag $asmbl_id\n" if ($DEBUG);
    }
    close (ATTFILE);

    foreach my $feat_id (keys %feat_hash) {
	if ($feat_id =~ /^CONTEXT[0-9]+:(.*)$/) {
	    my $real_feat_id = $1;
	    if (!defined $feat_hash{$real_feat_id}) {
		print STDERR "ERROR: $real_feat_id did not appear in the gene attribute file but $feat_id did!\n";
		$failed = 1;
	    }
	} elsif (!defined $feat_hash{$feat_id}->{'end5'}) {
	    print STDERR "ERROR: $feat_id appears in the protein file ($pep_file) but not in the gene attribute file ($att_file)!\n";
	    $failed = 1;
	}
    }
    if ($failed) {
	exit(1);
    }

    foreach $tag (keys %genome_hash) {
	foreach $asmbl_id (keys %{ $genome_hash{$tag} }) {
	    my $sort_by_end5_end3 = sub {
		my $end5a = $feat_hash{$a}->{'end5'};
		my $end5b = $feat_hash{$b}->{'end5'};
		my $end3a = $feat_hash{$a}->{'end3'};
		my $end3b = $feat_hash{$b}->{'end3'};
		my $mida = ($end5a + $end3a) / 2;
		my $midb = ($end5b + $end3b) / 2;
		
		if ($mida <=> $midb) {
		    return ($mida <=> $midb);
		} elsif ($end5a <=> $end5b) {
		    return ($end5a <=> $end5b);
		} elsif ($end3a <=> $end3b) {
		    return ($end3a <=> $end3b);
		} else {
		    die ("SORRY $a and $b have the same end5 ($end5a) and end3 ($end3a), please correct the gene_att file!\n");
		}
	    };
	    
	    @{ $genome_hash{$tag}{$asmbl_id} } = sort $sort_by_end5_end3 ( @{ $genome_hash{$tag}{$asmbl_id} } );
	    for my $index (0 .. $#{ $genome_hash{$tag}{$asmbl_id} }) {
		my $feat_name = $genome_hash{$tag}{$asmbl_id}->[$index];
		if ($feat_name =~ /^CONTEXT[0-9]+:(.*)$/) {
		    $genome_hash{$tag}{$asmbl_id}->[$index] = $1;
		    $genome_hash_context{$tag}{$asmbl_id}->{$index} = 1;
		    delete $feat_hash{$feat_name};
		    print STDERR "CONTEXT: $tag $asmbl_id $index = $feat_name ($1)\n" if ($DEBUG);
		} else {
		    $TagByPointer{$feat_name} = $index;
		    print STDERR "$tag $asmbl_id $index = $feat_name\n" if ($DEBUG);
		}
	    }
	}
    }
    return;
}

sub rebuild_genome_hash {
    #rebuild %genome_hash, %genome_hash_context, and %TagByPointer after deleting any features/genes/proteins
    foreach my $tag (keys %genome_hash) {
	foreach my $asmbl_id (keys %{ $genome_hash{$tag} }) {
	    my $cur_index = 0;
	    my $skip_index = 0;
	    foreach my $feat_id (@{ $genome_hash{$tag}{$asmbl_id} }) {
		if (!defined $feat_hash{$feat_id}) {
		    $skip_index++;
		} else {
		    if ($cur_index < $skip_index) {
			if ((defined $genome_hash_context{$tag}{$asmbl_id}->{$cur_index}) && (!defined $genome_hash_context{$tag}{$asmbl_id}->{$skip_index})) {
			    delete $genome_hash_context{$tag}{$asmbl_id}->{$cur_index};
			}
			if ((!defined $genome_hash_context{$tag}{$asmbl_id}->{$cur_index}) && (defined $genome_hash_context{$tag}{$asmbl_id}->{$skip_index})) {
			    $genome_hash_context{$tag}{$asmbl_id}->{$cur_index} = 1;
			}
			$genome_hash{$tag}{$asmbl_id}->[$cur_index] = $genome_hash{$tag}{$asmbl_id}->[$skip_index];
			if (!defined $genome_hash_context{$tag}{$asmbl_id}->{$skip_index}) {
			    $TagByPointer{$feat_id} = $cur_index;
			    print STDERR "$tag $asmbl_id $feat_id $cur_index $skip_index\n" if ($DEBUG);
			}
		    }
		    $cur_index++;
		    $skip_index++;
		}
	    }
	    $cur_index--;
	    $#{ $genome_hash{$tag}{$asmbl_id} } = $cur_index;
	    $cur_index++;
	    while ($cur_index < $skip_index) {
		if (defined $genome_hash_context{$tag}{$asmbl_id}->{$cur_index}) {
		    delete $genome_hash_context{$tag}{$asmbl_id}->{$cur_index};
		}
		$cur_index++;
	    }
	}
    }
    return;
}

sub select_max_scores_from_btab { # get tab-delimited BLAST results in either WUBLAST or NCBI -m8/-m9 formats
## new code - test for btab format (WUbtab or ncbi m8/m9)

    my @btab_line = (); # array variable to store split btab lines
    my $qid = ""; # query id
    my $sid = ""; # subject id (from database)
    my $qtag; # query genome tag
    my $stag; # subject genome tag
    my $score = ""; # BLAST bit score

    open (my $infile, "<", "$btabpath/$btabfile") || die ("ERROR: can't open file $btabpath/$btabfile\n");

    # this while look is interrogate the file to determine which tabular output style it is
    while (<$infile>)  {
	chomp;
	if (!/^#/) { # don't look at REMs
	    @btab_line = split(/\t/);
	    last; # we have out info, break the loop
	}
    }
    close($infile);

    if ($#btab_line >= "19") { # adjusted because Badgers btab.pl script only goes to e-value or perl column 19
	print STDERR "Detected WUBLAST-style btab file ...\n";
	$btab_style = 1; # WU
    }
    elsif ($#btab_line == "11") {
	print STDERR "NCBI BLAST (-m 8 or -m 9) option btab file detected ...\n";
	$btab_style = 0; # NCBI
    }
    else {
	die ("ERROR:  BLAST data must be either WUBLAST btab or NCBI BLAST -m8 or -m9 formats.\n");
    }
    ### process BLAST results ###
    open ($infile, "<", "$btabpath/$btabfile");
    while (<$infile>)  {
	chomp;
	@btab_line = split(/\t/);
	if ($btab_style) { # WU
	    $score = $btab_line[12];
	    $qid = $btab_line[0];
	    $sid = $btab_line[5];
	} else { # NCBI
	    if ($btab_line[0] =~ /^#/) { next;} # skip the lines beginning with #
	    $score = $btab_line[11];
	    $qid = $btab_line[0];
	    $sid = $btab_line[1];
	}
	if (!defined($qid) || !defined $feat_hash{$qid}) {print STDERR "WARNING!!! $qid is a feature identifier in the btab file but not in the gene attribute file skipping this entry!\n" if ($DEBUG); next;}
        if (!defined($sid) || !defined $feat_hash{$sid}) {print STDERR "WARNING!!! $sid is a feature identifier in the btab file but not in the gene attribute file skipping this entry!\n" if ($DEBUG); next;}
	$qtag = $FeatnameLookupTag_hash{$qid};
	$stag = $FeatnameLookupTag_hash{$sid};
	if (!defined($qtag) || !defined $TagIndex{$qtag}) {print STDERR "WARNING!!! $qtag is a genome tag in the btab file but not in the genome tag file skipping this entry!\n" if ($DEBUG); next;}
	if (!defined($stag) || !defined $TagIndex{$stag}) {print STDERR "WARNING!!! $stag is a genome tag in the btab file but not in the genome tag file skipping this entry!\n" if ($DEBUG); next;}
	print STDERR "$qtag $qid $stag $sid $score\n" if ($DEBUG);
	if (!defined $BestTagMatch{$qid}{$stag}) { # this assumes that the best Blast match for a query in a given subject genome appears first in the file
	    $BestTagMatch{$qid}->{$stag} = $score;
	} else {
	    if (!defined $SecondBestTagMatch{$qid}{$stag}) { # this assumes that the second best Blast match for a query in a given subject genome appears second in the file
		$SecondBestTagMatch{$qid}->{$stag} = $score;
	    }
	    if ($score > $BestTagMatch{$qid}->{$stag}) { # correct Best scores if sort assumption is violated
		$SecondBestTagMatch{$qid}->{$stag} = $BestTagMatch{$qid}->{$stag};
		$BestTagMatch{$qid}->{$stag} = $score;
	    } elsif ($score > $SecondBestTagMatch{$qid}->{$stag}) {
		$SecondBestTagMatch{$qid}->{$stag} =  $score;
	    }
	}
	#adding in recipricol direction in case blast results lack symmetry as they sometimes do
	#this will violate the sort assumptions but we can correct for that
	#this may mask problems with the all against all having missing data
	if (!defined $BestTagMatch{$sid}{$qtag}) {
	    $BestTagMatch{$sid}->{$qtag} = $score;
	} else {
	    if (!defined $SecondBestTagMatch{$sid}{$qtag}) {
		$SecondBestTagMatch{$sid}->{$qtag} = $score;
	    }
	    if ($score > $BestTagMatch{$sid}->{$qtag}) { # correct Best scores if sort assumption is violated
		$SecondBestTagMatch{$sid}->{$qtag} = $BestTagMatch{$sid}->{$qtag};
		$BestTagMatch{$sid}->{$qtag} = $score;
	    } elsif ($score > $SecondBestTagMatch{$sid}->{$qtag}) {
		$SecondBestTagMatch{$sid}->{$qtag} =  $score;
	    }
	}
    }
    close ($infile);
}

sub select_data_from_btab { # get tab-delimited BLAST results in either WUBLAST or NCBI -m8/-m9 formats
## new code - test for btab format (WUbtab or ncbi m8/m9)

    my @btab_line = (); # array variable to store split btab lines
    my $qmatch_length = ""; # stores the query length
    my $smatch_length = ""; # stores the subject length
    my $qid = ""; # query id
    my $sid = ""; # subject id (from database)
    my $qtag; # query genome tag
    my $stag; # subject genome tag
    my $qbegin = ""; # start query
    my $qend = ""; # end query
    my $sbegin = ""; # start subject
    my $send = ""; # end subject
    my $pid = ""; # percent identity
    my $evlu = ""; # e-value
    my $score = ""; # BLAST bit score
    my $qlength = ""; # size of query protein sequence
    my $slength = ""; # size of subject (database match) protein sequence

    open (my $infile, "<", "$btabpath/$btabfile") || die ("ERROR: can't open file $btabpath/$btabfile\n");

    ### process BLAST results ###
    open ($infile, "<", "$btabpath/$btabfile");
    while (<$infile>)  {
	chomp;
	@btab_line = split(/\t/);
	# same variables for both btab styles
	$qbegin = $btab_line[6];
	$qend = $btab_line[7];
	$sbegin = $btab_line[8];
	$send = $btab_line[9];
	
	if ($btab_style) { # WU
		
	    # adjusted because Badgers WU btab.pl script only goes to e-value or perl column 19
	    # ========================================================
	    # btab output for WUBLAST output
	    # column number Description (for Perl), add 1 for Unix
	    # 0       Query Sequence Name
	    # 1       Date of the Analysis
	    # 2       Query Sequence Length
	    # 3       Search Method  --  Blast family application name
	    # 4       Database Name
	    # 5       Subject Sequence Name  --  Database entry name
	    # 6       Start of alignment on query (5' nucleotide match in query)
	    # 7       End of alignment on query (3' nucleotide match in query)
	    # 8       Start of alignment on subject (5' nucleotide match in db hit)
	    # 9       End of alignment on subject (3' nucleotide match in db hit)
	    # 10      % Identity 
	    # 11      % Similarity 
	    # 12      Score (bits)
	    # 13      File Offset for Beginning of Alignment
	    # 14      File Offset for End of Alignment
	    # 15      Description (annotatioon)
	    # 16      Frame  --  1 through 6, or NULL
	    # 17      Query Strand  --  Plus, Minus or NULL
	    # 18      DB sequence length
	    # 19      Expect -- expected value
	    # 20      P-Value  --  Poisson ratio
	    # ========================================================
		
	    # propigate variables for WUBLAST
	    $qid = $btab_line[0];
	    $sid = $btab_line[5];
	    $pid = $btab_line[10];
	    $evlu = $btab_line[19];
	    $score = $btab_line[12];
	} else { # NCBI
	    # ========================================================
	    # btab output from NCBI blastn (-m 8) option:
	    # column number Description (for Perl), add 1 for Unix
	    # 0      Query_id
	    # 1	     subject_id (Hit from db)
	    # 2	     % Identity
	    # 3	     length of alignment
	    # 4	     number or mismatches
	    # 5	     number of gaps
	    # 6	     start of alignment on query (5' nucleotide match in query)
	    # 7	     end of alignment on query (3' nucleotide match in query)
	    # 8	     start of alignment on subject (5' nucleotide match in db hit)
	    # 9	     end of alignment on subject (3' nucleotide match in db hit)
	    # 10     e-value
	    # 11     score (bits)
	    # ========================================================
	    
	    # propigate variables for NCBI m8/m9
	    if ($btab_line[0] =~ /^#/) { next;} # skip the lines beginning with #
	    $qid = $btab_line[0];
	    $sid = $btab_line[1];
	    $pid = $btab_line[2];
	    $evlu = $btab_line[10];
	    $score = $btab_line[11];
	}
	### Generic processing ###
	print STDERR "$qid $sid $score\n" if ($DEBUG);
	if (!defined($qid) || !defined $feat_hash{$qid}) {print STDERR "WARNING!!! $qid is a feature identifier in the btab file but not in the gene attribute file skipping this entry!\n" if ($DEBUG); next;}
        if (!defined($sid) || !defined $feat_hash{$sid}) {print STDERR "WARNING!!! $sid is a feature identifier in the btab file but not in the gene attribute file skipping this entry!\n" if ($DEBUG); next;}
	$qtag = $FeatnameLookupTag_hash{$qid};
	$stag = $FeatnameLookupTag_hash{$sid};
	if (!defined($qtag) || !defined $TagIndex{$qtag}) {print STDERR "WARNING!!! $qtag is a genome tag in the btab file but not in the genome tag file skipping this entry!\n" if ($DEBUG); next;}
	if (!defined($stag) || !defined $TagIndex{$stag}) {print STDERR "WARNING!!! $stag is a genome tag in the btab file but not in the genome tag file skipping this entry!\n" if ($DEBUG); next;}
	$qlength = $feat_hash{$qid}->{'length'};
	$slength = $feat_hash{$sid}->{'length'};
	$qmatch_length = ((abs($qbegin - $qend) + 1)/$qlength)*100; 
	$smatch_length = ((abs($sbegin - $send) + 1)/$slength)*100;
	if (!$feat_hash{$qid}->{'bit'})  {
	    $feat_hash{$qid}->{'bit'} = 1; # used to determine number of orfs processed
	    $orf_counter{$qtag}->{'raw'}++;  # increment the total orf counter (we will get the total # of orfs this way)
	}
	if (!defined $relationship_hash{$qid}) { # this assumes that the best Blast match (presumably to itself) for a query appears first in the file
	    $query_evalue_cutoff{$qid} = $evlu * 10.0;
	}
	if ($pid < $percentid)  {
	    print STDERR "Skipping %id $pid < $percentid Query: $qid X Subject $sid\n" if ($DEBUG);
	    next;
	}
	if (($evlu > $evalue) && ($evlu > $query_evalue_cutoff{$qid}))  {
	    print STDERR "Skipping $evlu > $evalue Query: $qid X Subject $sid = $pid\n" if ($DEBUG);
	    next;
	}
	if (($qmatch_length < $min_hit_length) || ($smatch_length < $min_hit_length))  {
	    print STDERR "Skipping match length $qmatch_length , $smatch_length < $min_hit_length Query: $qid X Subject $sid = $pid\n" if ($DEBUG);
	    next;
	}
	#if (($qmatch_length < 50.0) && ($smatch_length < 50.0))  { # the match length must be greater than 1/2 the length of the shorter protein to ignore domain matches
	#print STDERR "Skipping match length $qmatch_length , $smatch_length both < 50.0 Query: $qid X Subject $sid = $pid\n" if ($DEBUG);
	#next;
	#} Didn't help - but did hurt finding some valid fragments
	if (!defined $BestTagMatch{$qid}{$stag}) {
	    print STDERR "ERROR: $qid and genome $stag had no score the first time through!\n";
	    exit(1);
	}
	if ((0.8 * $BestTagMatch{$qid}->{$stag}) > $score) {
	    print STDERR "Skipping score $score < 80% of max $BestTagMatch{$qid}->{$stag} Query: $qid X Subject $sid = $pid\n" if ($DEBUG);
	    next; #do not keep matches scoring less than 80% of the best score
	}
	if (defined $relationship_hash{$qid} && defined $relationship_hash{$qid}{$sid}) {
	    if ($relationship_hash{$qid}{$sid}->{'score'} >= $score) { # ignore lower scores for the same two proteins
		next;
	    }
	} else {
	    $Qbytaghash{$qid}{$stag}->{$sid} = $relationship_hash{$qid}{$sid} = {}; #have Qbytaghash and relationship_hash reference the same hash
	}
	    # make sure $qbegin is smaller than $qend and $sbegin is smaller than $send
	my $temp_val;
	if ($qbegin > $qend) { #swap
	    $temp_val = $qend;
	    $qend = $qbegin;
	    $qbegin = $temp_val;
	}
	if ($sbegin > $send) { #swap
	    $temp_val = $send;
	    $send = $sbegin;
	    $sbegin = $temp_val;
	}
	$relationship_hash{$qid}{$sid}->{'id'} = $pid;
	$relationship_hash{$qid}{$sid}->{'eval'} = $evlu;
	$relationship_hash{$qid}{$sid}->{'score'} = $score;
	$relationship_hash{$qid}{$sid}->{'min_query'} = $qbegin;
	$relationship_hash{$qid}{$sid}->{'max_query'} = $qend;
	$relationship_hash{$qid}{$sid}->{'min_sub'} = $sbegin;
	$relationship_hash{$qid}{$sid}->{'max_sub'} = $send;
	$relationship_hash{$qid}{$sid}->{'best'} = 0;
	$relationship_hash{$qid}{$sid}->{'bibest'} = 0;
	$relationship_hash{$qid}{$sid}->{'synbest'} = 0;
	$relationship_hash{$qid}{$sid}->{'synbibest'} = 0;
	$relationship_hash{$qid}{$sid}->{'CGN_bibest'} = 0;
	$relationship_hash{$qid}{$sid}->{'full'} = 0;
	$relationship_hash{$qid}{$sid}->{'anchor'} = 0;
	$relationship_hash{$qid}{$sid}->{'extend'} = 0;
	print STDERR "Query: $qid X Subject $sid = $relationship_hash{$qid}{$sid}->{'id'}\n" if ($DEBUG);
	#determine full length indicating blast matches
	if ($pid >= $fspercentid) { # only use matches above the defined percent identity cutoff
	    if (($qbegin <= $max_missing_aa + 1) && ($qend + $max_missing_aa >= $feat_hash{$qid}->{'length'}) && ($qmatch_length >= 80)) { # covers full length of query
		if (($sbegin <= $max_missing_aa + 1) && ($send + $max_missing_aa >= $feat_hash{$sid}->{'length'}) && ($smatch_length >= 80)) { #covers full length of subject
		    $feat_hash{$qid}->{'full'}++; #not a protein fragment if almost a full length match
		    $feat_hash{$sid}->{'full'}++; #this is also true for the subject but we may double count
		    $relationship_hash{$qid}{$sid}->{'full'} = 1;
		}
	    }
	}
	if (!defined $Tagbyfeatnamehash{$qtag}{$qid})  {
	    $Tagbyfeatnamehash{$qtag}{$qid} = 1;
	    $orf_counter{$qtag}->{'used'}++;
	} elsif ($Tagbyfeatnamehash{$qtag}{$qid} == 0) { #if qid had only been seen as a target (sid) before change to seen as a query (qid)
		$Tagbyfeatnamehash{$qtag}{$qid} = 1;
	}
	if (!defined $Tagbyfeatnamehash{$stag}{$sid})  {
	    $Tagbyfeatnamehash{$stag}{$sid} = 0;
	    $orf_counter{$stag}->{'used'}++;
	}
    }
    close ($infile);

    open (OUTMISSING, ">$outprefix" . "missing_blast_results.txt");
    my $changed = 0;
    foreach my $qid (keys %feat_hash)  { # go through all featnames
	my $qtag = $FeatnameLookupTag_hash{$qid};
	if (!defined $Tagbyfeatnamehash{$qtag}{$qid})  {#output feat_names in the attribute file but missing from the blast results, and cleanup hashes to ignore these proteins
	    $changed = 1;
	    print OUTMISSING "$qtag:$qid\n";
	    delete $feat_hash{$qid};
	    delete $FeatnameLookupTag_hash{$qid};
	    delete $AssemblyLookup_hash{$qid};
	    delete $TagByPointer{$qid};
	} elsif ($Tagbyfeatnamehash{$qtag}{$qid} == 0) { #if qid had only been seen as a target
	    print OUTMISSING "$qtag:$qid only search result\n";
	}
    }
    if ($changed) {
	&rebuild_genome_hash; #need to remove the deleted features/genes/proteins
    }
    close (OUTMISSING);
    print STDERR "check that TagByPointer is consistent\n" if ($DEBUG);
#check that TagByPointer is consistent
    foreach $qtag (@tag_array)  {  # start looping through genomes by order in tag file
	foreach $qid (keys %{ $Tagbyfeatnamehash{$qtag} } ) { # go through featnames of each genome to look for matches in other genomes
	    foreach $stag (keys %{ $Qbytaghash{$qid} } ) {# if query protein matches anything in subject genome, lets drill through each match
		foreach $sid (keys %{ $Qbytaghash{$qid}{$stag} } ) {
		    my $qasmbl_id = $AssemblyLookup_hash{$qid};
		    if (!defined $qasmbl_id) {
			die ("ERROR: asmbl_id undefined for $qtag $qid ($stag $sid)\n");
		    }
		    my $sasmbl_id = $AssemblyLookup_hash{$sid};
		    if (!defined $sasmbl_id) {
			die ("ERROR: asmbl_id undefined for ($qtag $qid) $stag $sid\n");
		    }
		    my $qindex = $TagByPointer{$qid};
		    if (!defined $qindex) {
			die ("ERROR: TagByPointer undefined for $qtag $qid ($stag $sid)\n");
		    }
		    my $sindex = $TagByPointer{$sid};
		    if (!defined $sindex) {
			die ("ERROR: TagByPointer undefined for ($qtag $qid) $stag $sid\n");
		    }
		    if ($genome_hash{$qtag}{$qasmbl_id}->[$qindex] ne $qid) {
			die ("ERROR: Inconsistent $genome_hash{$qtag}{$qasmbl_id}->[$qindex] $qid for $qtag $qid ($stag $sid)\n");
		    }
		    if ($genome_hash{$stag}{$sasmbl_id}->[$sindex] ne $sid) {
			die ("ERROR: Inconsistent $genome_hash{$stag}{$sasmbl_id}->[$sindex] $sid for ($qtag $qid) $stag $sid\n");
		    }
		}
	    }
	}
    }


#check that BestMatchTag and SecondBestMatchTag are consistent
#also determine if self-hit exists and create self-hit record if not, no score should be above the self-hit so reset them
    my $cum_num_qids = 0; #calculate average number of proteins in a genome
    foreach my $qtag (@tag_array)  {  # start looping through genomes by order in tag file
	my $num_qids = 0;
	my $num_no_selfhits = 0;
	foreach my $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome to look for matches in other genomes
	    $num_qids++;
	    my $max_score = 0;
	    my $selfhit = 1;
	    my $selfhit_score;
	    if (defined $relationship_hash{$qid} && defined $relationship_hash{$qid}{$qid}) {
		$selfhit_score = $relationship_hash{$qid}{$qid}->{'score'};
		if (!defined $BestTagMatch{$qid}{$qtag}) {
		    die ("ERROR: selfhit was found previously but not here!\n");
		}
	    } else {
		print STDERR "WARNING: no self hit Blast score for $qtag : $qid\n" if ($DEBUG);
		$selfhit = 0;
		$num_no_selfhits++;
	    }
	    if (!defined $BestTagMatch{$qid}) {
		die ("ERROR: no  Blast scores for $qtag : $qid when some were found in previous loop\n");
	    }
	    # Find maximum blast score to use for selfhit score if there is no selfhit score
	    if (!$selfhit) {
		foreach my $stag (@tag_array) { # loop through all genomes
		    if (defined $BestTagMatch{$qid}{$stag} && ($BestTagMatch{$qid}{$stag} > $max_score)) {
			$max_score = $BestTagMatch{$qid}{$stag};
		    }
		}
		$selfhit_score = $max_score;
	    }
	    foreach my $stag (@tag_array) { # loop through all genomes
		if (defined $SecondBestTagMatch{$qid}{$stag} && !defined $BestTagMatch{$qid}{$stag}) {
		    die ("ERROR: Second Best Blast score for $qtag : $qid to genome $stag exists but not Best!\n");
		}
		if (!defined $BestTagMatch{$qid}{$stag}) {
		    next;
		}
		if ($BestTagMatch{$qid}{$stag} > $selfhit_score) {
		    $BestTagMatch{$qid}{$stag} = $selfhit_score;
		    if (defined $SecondBestTagMatch{$qid}{$stag} && ($SecondBestTagMatch{$qid}{$stag} > $selfhit_score)) {
			$SecondBestTagMatch{$qid}{$stag} = $selfhit_score;
		    }
		}
	    }
	    if (!$selfhit) {
		$BestTagMatch{$qid}{$qtag} = $max_score;
		if ((defined $Qbytaghash{$qid}{$qtag}->{$qid}) || (defined $relationship_hash{$qid}{$qid})) {
		    die ("ERROR: selfhit defined here ($qtag $qid) but not above\n");
		}
		$Qbytaghash{$qid}{$qtag}->{$qid} = $relationship_hash{$qid}{$qid} = {}; #have Qbytaghash and relationship_hash reference the same hash
		$relationship_hash{$qid}{$qid}->{'score'} = $max_score;
		$relationship_hash{$qid}{$qid}->{'id'} = 100;
		$relationship_hash{$qid}{$qid}->{'eval'} = 0;
		$relationship_hash{$qid}{$qid}->{'min_query'} = 1;
		$relationship_hash{$qid}{$qid}->{'max_query'} = $feat_hash{$qid}->{'length'};
		$relationship_hash{$qid}{$qid}->{'min_sub'} = 1;
		$relationship_hash{$qid}{$qid}->{'max_sub'} = $feat_hash{$qid}->{'length'};
		$relationship_hash{$qid}{$qid}->{'best'} = 0;
		$relationship_hash{$qid}{$qid}->{'bibest'} = 0;
		$relationship_hash{$qid}{$qid}->{'synbest'} = 0;
		$relationship_hash{$qid}{$qid}->{'synbibest'} = 0;
		$relationship_hash{$qid}{$qid}->{'CGN_bibest'} = 0;
		$relationship_hash{$qid}{$qid}->{'full'} = 0;
		$relationship_hash{$qid}{$qid}->{'anchor'} = 0;
		$relationship_hash{$qid}{$qid}->{'extend'} = 0;
	    }
	}
	$cum_num_qids += $num_qids;
	if ($num_no_selfhits > 0) {
	    print STDERR "WARNING: Number of missing Blast selfhits for genome $qtag is $num_no_selfhits\n";
	}
    }
    my $ave_num_qids = 0;
    if ($genome_number > 0) {
	$ave_num_qids = $cum_num_qids / $genome_number;
    }

    print STDERR "reduce any scores above self score to self score\n" if ($DEBUG);
#reduce any scores above self score to self score
    foreach my $qid (keys %feat_hash)  { # go through all featnames
	if (!defined $relationship_hash{$qid}) {
	    print STDERR "WARNING: no matches to any other gene for $qid!\n";
	    next;
	}
	if (!defined $relationship_hash{$qid}{$qid}) {
	    die ("ERROR: $qid is missing a selfhit after forcing one!\n");
	}
	foreach my $sid (keys %{ $relationship_hash{$qid} })  { # go through all featnames
	    if ($relationship_hash{$qid}{$sid}->{'score'} > $relationship_hash{$qid}{$qid}->{'score'}) {
		print STDERR "WARNING ($qid $sid $relationship_hash{$qid}{$sid}->{'score'}) > selfhit $relationship_hash{$qid}{$qid}->{'score'}\n" if ($DEBUG);
		$relationship_hash{$qid}{$sid}->{'score'} = $relationship_hash{$qid}{$qid}->{'score'};
	    }
	}
    }

#check that genomes have a reasonable amount of proteins
    foreach my $qtag (@tag_array)  {  # start looping through genomes by order in tag file
	my $num_qids = 0;
	foreach my $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome
	    $num_qids++;
	}
	print STDERR "genome: $qtag has $num_qids valid genes";
	if ($num_qids < (0.75 * $ave_num_qids)) {
	    print STDERR " WARNING!!! This is much lower than the average number of genes $ave_num_qids";
	}
	print STDERR "\n";
    }
}

sub calc_BSR  { # subroutine to generate the average (both directions) Blast Score Ratio (BSR) and populate Qbytaghash and relationship_hash with BSR

    my $qid = ""; # query feat_name
    my $sid = ""; # subject feat_name
    my $qtag = ""; # query genome tag
    my $stag = ""; # subject tag

    foreach $qtag (@tag_array)  {  # start looping through by order in tag file (reference is first)
	foreach $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome to look for matches in other genomes
	    foreach $stag (@tag_array)  {
		if (defined $Qbytaghash{$qid}{$stag}) { # if query protein matches anything in subject genome, lets drill through each match
		    foreach $sid (keys %{ $Qbytaghash{$qid}{$stag} })  { # need to check value of $sid here (see comment in previous subroutine)
			# calculate the Blast Score Ratio (BSR) of each relationship
			if (!defined $relationship_hash{$sid}{$qid}) {
			    $Qbytaghash{$sid}{$qtag}->{$qid} = $relationship_hash{$sid}{$qid} = {}; #have Qbytaghash and relationship_hash reference the same hash
			    #setting these values to make a symmetric entry when the match in one direction was discarded
			    if ($relationship_hash{$qid}{$sid}->{'score'} > $relationship_hash{$sid}{$sid}->{'score'}) { # prevents norm_score being greater than 1
				$relationship_hash{$sid}{$qid}->{'score'} = $relationship_hash{$sid}{$sid}->{'score'};
			    } else {
				$relationship_hash{$sid}{$qid}->{'score'} = $relationship_hash{$qid}{$sid}->{'score'};
			    }
			    $relationship_hash{$sid}{$qid}->{'id'} = $relationship_hash{$qid}{$sid}->{'id'};
			    $relationship_hash{$sid}{$qid}->{'eval'} = $relationship_hash{$qid}{$sid}->{'eval'};
			    $relationship_hash{$sid}{$qid}->{'min_query'} = $relationship_hash{$qid}{$sid}->{'min_sub'};
			    $relationship_hash{$sid}{$qid}->{'max_query'} = $relationship_hash{$qid}{$sid}->{'max_sub'};
			    $relationship_hash{$sid}{$qid}->{'min_sub'} = $relationship_hash{$qid}{$sid}->{'min_query'};
			    $relationship_hash{$sid}{$qid}->{'max_sub'} = $relationship_hash{$qid}{$sid}->{'max_query'};
			    $relationship_hash{$sid}{$qid}->{'best'} = 0;
			    $relationship_hash{$sid}{$qid}->{'bibest'} = 0;
			    $relationship_hash{$sid}{$qid}->{'synbest'} = 0;
			    $relationship_hash{$sid}{$qid}->{'synbibest'} = 0;
			    $relationship_hash{$sid}{$qid}->{'CGN_bibest'} = 0;
			    $relationship_hash{$sid}{$qid}->{'anchor'} = 0;
			    $relationship_hash{$sid}{$qid}->{'extend'} = 0;
			    $relationship_hash{$sid}{$qid}->{'full'} = $relationship_hash{$qid}{$sid}->{'full'};
			}
			if (!defined $relationship_hash{$qid}{$qid}) {
			    die ("ERROR: somehow we got to calc_BSR without a selfhit for $qid!\n");
			}
			my $norm_score;
			$relationship_hash{$qid}{$sid}->{'norm_score'} = $norm_score = $relationship_hash{$qid}{$sid}->{'score'} / $relationship_hash{$qid}{$qid}->{'score'};
			if (($norm_score > 1) || ($norm_score < 0)) {
			    die ("ERROR: normalized BSR score for $qtag:$qid,$stag:$sid is not between 0-1 ($norm_score)\n");
			}
			if (!defined $relationship_hash{$sid}{$sid}) {
			    die ("ERROR: somehow we got to calc_BSR without a selfhit for $sid!\n");
			}
			$relationship_hash{$sid}{$qid}->{'norm_score'} = $norm_score = $relationship_hash{$sid}{$qid}->{'score'} / $relationship_hash{$sid}{$sid}->{'score'};
			if (($norm_score > 1) || ($norm_score < 0)) {
			    die ("ERROR: normalized BSR score for $stag:$sid;$qtag:$qid is not between 0-1 ($norm_score)\n");
			}
			#BSR (Blast Score Ratio) is mean of two normalized scores
			$Qbytaghash{$qid}{$stag}{$sid}->{'BSR'} = $Qbytaghash{$sid}{$qtag}{$qid}->{'BSR'} = ($relationship_hash{$qid}{$sid}->{'norm_score'} + $relationship_hash{$sid}{$qid}->{'norm_score'})/2 ; # calc BSR GGS need to keep $Qbytaghash symmetrical for frameshifts
			#alternate choice of BSR to be max of two normalized scores - this did not work but we do use the max of the norm scores in a few places
			#if ($relationship_hash{$qid}{$sid}->{'norm_score'} > $relationship_hash{$sid}{$qid}->{'norm_score'}) {
			#    $Qbytaghash{$qid}{$stag}->{$sid}->{'BSR'} = $Qbytaghash{$sid}{$qtag}{$qid}->{'BSR'} = $relationship_hash{$qid}{$sid}->{'norm_score'}; # calc BSR GGS need to keep $Qbytaghash symmetrical for frameshifts
			#} else {
			#    $Qbytaghash{$qid}{$stag}->{$sid}->{'BSR'} = $Qbytaghash{$sid}{$qtag}{$qid}->{'BSR'} = $relationship_hash{$sid}{$qid}->{'norm_score'}; # calc BSR GGS need to keep $Qbytaghash symmetrical for frameshifts
			#}
		    }
		}
	    }
	}
    }
			
    foreach my $qid (keys %relationship_hash)  { # go through all featnames
	foreach my $sid (keys %{ $relationship_hash{$qid} } )  {
	    if ((!defined $relationship_hash{$qid}{$sid}->{'BSR'}) || (!defined $relationship_hash{$sid}{$qid}->{'BSR'})) {
		print STDERR "ERROR!!! BSR undefined for $qid $sid in relationship_hash\n";
	    }
	}
    }
}

sub calc_bibest  { # subroutine to find bidrectional best blast matches and populate Qbytaghash and relationship_hash with them

    my $qid = ""; # query feat_name
    my $sid = ""; # subject feat_name
    my $qtag = ""; # query genome tag
    my $stag = ""; # subject tag

    #mark best blast hit per genome per qid
    foreach $qtag (@tag_array)  {  # start looping through by order in tag file (reference is first)
	foreach $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome to look for matches in other genomes
	    if (!defined $Qbytaghash{$qid}) {
		next;# need to do this to prevent keys %{ $Qbytaghash{$qid}{$stag} } causing $Qbytaghash{$qid}{$stag} to be defined later
	    }
	    foreach $stag (keys %{ $Qbytaghash{$qid} })  {
		if ($qtag eq $stag) {
		    next; #skip matches within the same genome
		}
		my $best_score = 0;
		my $second_score = 0;
		my $found_best = 0;
		my $best_sid;
		foreach $sid (sort {$Qbytaghash{$qid}{$stag}{$b}->{'BSR'} <=> $Qbytaghash{$qid}{$stag}{$a}->{'BSR'}} keys %{ $Qbytaghash{$qid}{$stag} })  {
		    if (!$found_best) {
			$best_score = $relationship_hash{$qid}{$sid}->{'BSR'};
			$found_best = 1;
			$best_sid = $sid;
			next;
		    }
		    $second_score = $relationship_hash{$qid}{$sid}->{'BSR'};
		    last; #quit after finding second best blast match
		}
		if ($found_best && ((0.95 * $best_score) > $second_score)) {
		    $relationship_hash{$qid}{$best_sid}->{'best'} = 1;
		}
	    }
	}
    }

    #mark bidirectionally best blast hit per genome per qid if it exists
    foreach $qid (keys %relationship_hash)  { # go through all featnames
	foreach $sid (keys %{ $relationship_hash{$qid} } )  {
		if ($relationship_hash{$qid}{$sid}->{'best'} && $relationship_hash{$sid}{$qid}->{'best'}) {
		    $relationship_hash{$qid}{$sid}->{'bibest'}  = 1;
		}
	}
    }
			
}

sub print_pairwise_matrix  { # subroutine to print similarity or distance matrices in a few formats

    my ($matrix, $outfile, $scale, $offset, $upper, $phylip_dist) = @_;

    if ($phylip_dist) {
	print $outfile "$genome_number\n";
    } else {
	foreach my $tag (@tag_array)  {
	    printf $outfile "\t%7s", substr ($tag, -7);
	}
	print $outfile "\n";
    }
    foreach my $row_index (0 .. $#tag_array)  {
	printf $outfile "%-7s", substr ($tag_array[$row_index], -7);
	my $allowed_per_line = 1;
	foreach my $col_index (0 .. $#tag_array)  {
	    if ($phylip_dist) {
		if ($allowed_per_line == 6) {
		    print $outfile "\n";
		    $allowed_per_line = 0;
		} else {
		    print $outfile " ";
		}
	    } else {
		print $outfile "\t";
	    }
	    $allowed_per_line++;
	    if (!$upper || ($col_index > $row_index)) {
		printf $outfile "%7.2f", ($offset + ($scale * $matrix->[$row_index][$col_index]));
	    }
	}
	print $outfile "\n";
    }

}

sub calc_score_histograms  { # subroutine to calculate and output score histograms and score cutoffs

    my @tmp_tag_array = @tag_array;
    my %self_hist = ();

    open (OUTHISTOGRAM, ">$outprefix" . "histograms.txt") if ($histogramfile);

    foreach my $qtag (@tag_array)  {  # start looping through genomes by order in tag file
	$self_hist{$qtag} = [];
	foreach my $index (0..100) { # initialize self histogram
	    $self_hist{$qtag}[$index] = 0;
	}
	my $num_qids = 0;
	foreach my $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome to look for matches to paralogs in the same genome
	    $num_qids++;
	    my $match_found = 0;
	    foreach my $sid (sort {$Qbytaghash{$qid}{$qtag}{$b}->{'norm_score'} <=> $Qbytaghash{$qid}{$qtag}{$a}->{'norm_score'}} keys %{ $Qbytaghash{$qid}{$qtag} })  {
		if ($qid eq $sid) { # skip self match
		    next;
		}
		$self_hist{$qtag}[int (100 * $relationship_hash{$qid}{$sid}->{'norm_score'})]++;
		$match_found = 1;
		last; #quit after finding best nonself match
	    }
	    if (!$match_found) {
		if (!defined $SecondBestTagMatch{$qid}{$qtag}) {
		    $self_hist{$qtag}[0]++;
		} else {
		$self_hist{$qtag}[int (100 * ($SecondBestTagMatch{$qid}{$qtag} / $BestTagMatch{$qid}{$qtag}))]++;
		}
	    }
	}
	if ($num_qids > 0) { #do not divide by zero when there are no matches
	    foreach my $index (0..100) { # normalize self histogram
		$self_hist{$qtag}[$index] /= $num_qids;
	    }
	}
	if ($histogramfile) {
	    print OUTHISTOGRAM "$qtag:$qtag";
	    foreach my $index (0..100) { # print self histogram
		print OUTHISTOGRAM "\t$self_hist{$qtag}[$index]";
	    }
	    print OUTHISTOGRAM "\n";
	}
    }
    #initialize global arrays we are about to compute
    foreach my $index1 (0 .. $#tag_array) {
	foreach my $index2 (0 .. $#tag_array) {
	    $paralog_cutoff[$index1][$index2] = 0;
	    $ave_per_id[$index1][$index2] = 0;
	    $mean_max_BSR[$index1][$index2] = 0;
	}
    }

    foreach my $qtag (@tag_array)  {  # start looping through genomes by order in tag file
	my $qTagIndex = $TagIndex{$qtag};
	shift(@tmp_tag_array);
	foreach my $stag (@tmp_tag_array) { # loop through all other genomes
	    my $sTagIndex = $TagIndex{$stag};
	    my @qtag_best_hist;
	    my @qtag_second_hist;
	    my @qtag_good_ortho_hist;
	    my @qtag_good_second_hist;
	    my @stag_best_hist;
	    my @stag_second_hist;
	    my @stag_good_ortho_hist;
	    my @stag_good_second_hist;
	    my $qcum_per_id = 0;
	    my $scum_per_id = 0;
	    my $qcum_max_BSR = 0;
	    my $scum_max_BSR = 0;
	    foreach my $index (0..100) { # initialize histograms
		$qtag_best_hist[$index] = 0;
		$qtag_second_hist[$index] = 0;
		$qtag_good_ortho_hist[$index] = 0;
		$qtag_good_second_hist[$index] = 0;
		$stag_best_hist[$index] = 0;
		$stag_second_hist[$index] = 0;
		$stag_good_ortho_hist[$index] = 0;
		$stag_good_second_hist[$index] = 0;
	    }
	    my $num_qids = 0;
	    my $good_qids = 0;
	    foreach my $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome to look for matches in other genomes
		$num_qids++;
		if (!defined $Qbytaghash{$qid}) {
		    $qtag_best_hist[0]++;
		    $qtag_second_hist[0]++;
		    $qtag_good_ortho_hist[0]++;
		    $qtag_good_second_hist[0]++;
		    next;# need to do this to prevent keys %{ $Qbytaghash{$qid}{$stag} } causing $Qbytaghash{$qid}{$stag} to be defined later
		}
		if (!defined $Qbytaghash{$qid}{$stag}) {
		    $qtag_best_hist[0]++;
		    $qtag_second_hist[0]++;
		    $qtag_good_ortho_hist[0]++;
		    $qtag_good_second_hist[0]++;
		    next;
		}
		my $match_found = 0;
		my $best_found = 0;
		my $good_found = 0;
		foreach my $sid (sort {$Qbytaghash{$qid}{$stag}{$b}->{'norm_score'} <=> $Qbytaghash{$qid}{$stag}{$a}->{'norm_score'}} keys %{ $Qbytaghash{$qid}{$stag} })  {
		    if (!$best_found) {
			my $hist_index = int (100 * $relationship_hash{$qid}{$sid}->{'norm_score'});
			$qtag_best_hist[$hist_index]++;
			$best_found = 1;
			if ($relationship_hash{$qid}{$sid}->{'CGN_bibest'} > $CGN_window_size) {
			    $qtag_good_ortho_hist[$hist_index]++;
			    $qcum_per_id += $relationship_hash{$qid}{$sid}->{'id'};
			    my $max_BSR = $relationship_hash{$qid}{$sid}->{'norm_score'};
			    if ($relationship_hash{$sid}{$qid}->{'norm_score'} > $max_BSR) {
				$max_BSR = $relationship_hash{$sid}{$qid}->{'norm_score'};
			    }
			    $qcum_max_BSR += $max_BSR;
			    $good_qids++;
			    $good_found = 1;
			} else {
			    $qtag_good_ortho_hist[0]++;
			}
			next;
		    }
		    $qtag_second_hist[int (100 * $relationship_hash{$qid}{$sid}->{'norm_score'})]++;
		    if ($good_found) {
			$qtag_good_second_hist[int (100 * $relationship_hash{$qid}{$sid}->{'norm_score'})]++;
		    } else {
			$qtag_good_second_hist[0]++;
		    }
		    $match_found = 1;
		    last; #quit after finding second best match
		}
		if (!$best_found) {
		    $qtag_best_hist[0]++;
		    $qtag_good_ortho_hist[0]++;
		}
		if (!$match_found) {
		    if (!defined $SecondBestTagMatch{$qid}{$stag}) {
			$qtag_second_hist[0]++;
			$qtag_good_second_hist[0]++;
		    } else {
			my $best_score;
			if (!defined $BestTagMatch{$qid}{$qtag}) {
			    $best_score = $BestTagMatch{$qid}{$stag};
			} else {
			    $best_score = $BestTagMatch{$qid}{$qtag};
			}
			$qtag_second_hist[int (100 * ($SecondBestTagMatch{$qid}{$stag} / $best_score))]++;
			if ($good_found) {
			    $qtag_good_second_hist[int (100 * ($SecondBestTagMatch{$qid}{$stag} / $best_score))]++;
			} else {
			    $qtag_good_second_hist[0]++;
			}
		    }
		}
	    }
	    if ($num_qids > 0) { #do not divide by zero when there are no matches
		foreach my $index (0..100) { # normalize histograms
		    $qtag_best_hist[$index] /= $num_qids;
		    $qtag_second_hist[$index] /= $num_qids;
		    $qtag_good_ortho_hist[$index] /= $num_qids;
		    $qtag_good_second_hist[$index] /= $num_qids;
		}
	    }

	    my $num_sids = 0;
	    my $good_sids = 0;
	    foreach my $sid ( keys %{ $Tagbyfeatnamehash{$stag} } )  { # go through featnames of query genome to look for matches in other genomes
		$num_sids++;
		if (!defined $Qbytaghash{$sid}) {
		    $stag_best_hist[0]++;
		    $stag_second_hist[0]++;
		    $stag_good_ortho_hist[0]++;
		    $stag_good_second_hist[0]++;
		    next;# need to do this to prevent keys %{ $Qbytaghash{$qid}{$stag} } causing $Qbytaghash{$qid}{$stag} to be defined later
		}
		if (!defined $Qbytaghash{$sid}{$qtag}) {
		    $stag_best_hist[0]++;
		    $stag_second_hist[0]++;
		    $stag_good_ortho_hist[0]++;
		    $stag_good_second_hist[0]++;
		    next;
		}
		my $match_found = 0;
		my $best_found = 0;
		my $good_found = 0;
		foreach my $qid (sort {$Qbytaghash{$sid}{$qtag}{$b}->{'norm_score'} <=> $Qbytaghash{$sid}{$qtag}{$a}->{'norm_score'}} keys %{ $Qbytaghash{$sid}{$qtag} })  {
		    if (!$best_found) {
			my $hist_index = int (100 * $relationship_hash{$sid}{$qid}->{'norm_score'});
			$stag_best_hist[$hist_index]++;
			$best_found = 1;
			if ($relationship_hash{$sid}{$qid}->{'CGN_bibest'} > $CGN_window_size) {
			    $stag_good_ortho_hist[$hist_index]++;
			    $scum_per_id += $relationship_hash{$sid}{$qid}->{'id'};
			    my $max_BSR = $relationship_hash{$sid}{$qid}->{'norm_score'};
			    if ($relationship_hash{$qid}{$sid}->{'norm_score'} > $max_BSR) {
				$max_BSR = $relationship_hash{$qid}{$sid}->{'norm_score'};
			    }
			    $scum_max_BSR += $max_BSR;
			    $good_sids++;
			    $good_found = 1;
			} else {
			    $stag_good_ortho_hist[0]++;
			}
			next;
		    }
		    $stag_second_hist[int (100 * $relationship_hash{$sid}{$qid}->{'norm_score'})]++;
		    if ($good_found) {
			$stag_good_second_hist[int (100 * $relationship_hash{$sid}{$qid}->{'norm_score'})]++;
		    } else {
			$stag_good_second_hist[0]++;
		    }
		    $match_found = 1;
		    last; #quit after finding second best match
		}
		if (!$best_found) {
		    $stag_best_hist[0]++;
		    $stag_good_ortho_hist[0]++;
		}
		if (!$match_found) {
		    if (!defined $SecondBestTagMatch{$sid}{$qtag}) {
			$stag_second_hist[0]++;
			$stag_good_second_hist[0]++;
		    } else {
			my $best_score;
			if (!defined $BestTagMatch{$sid}{$stag}) {
			    $best_score = $BestTagMatch{$sid}{$qtag};
			} else {
			    $best_score = $BestTagMatch{$sid}{$stag};
			}
			$stag_second_hist[int (100 * ($SecondBestTagMatch{$sid}{$qtag} / $best_score))]++;
			if ($good_found) {
			    $stag_good_second_hist[int (100 * ($SecondBestTagMatch{$sid}{$qtag} / $best_score))]++;
			} else {
			    $stag_good_second_hist[0]++;
			}
		    }
		}
	    }
	    if ($num_sids > 0) { #do not divide by zero when there are no matches
		foreach my $index (0..100) { # normalize histograms
		    $stag_best_hist[$index] /= $num_sids;
		    $stag_second_hist[$index] /= $num_sids;
		    $stag_good_ortho_hist[$index] /= $num_sids;
		    $stag_good_second_hist[$index] /= $num_sids;
		}
	    }

	    if ($histogramfile) {
		print OUTHISTOGRAM "$qtag:$stag:good_CGN";
		foreach my $index (0..100) { # print qtag good ortho histogram
		    print OUTHISTOGRAM "\t$qtag_good_ortho_hist[$index]";
		}
		print OUTHISTOGRAM "\n";
		print OUTHISTOGRAM "$qtag:$stag:second_CGN";
		foreach my $index (0..100) { # print qtag good ortho histogram
		    print OUTHISTOGRAM "\t$qtag_good_second_hist[$index]";
		}
		print OUTHISTOGRAM "\n";
		print OUTHISTOGRAM "$qtag:$stag:best";
		foreach my $index (0..100) { # print qtag best histogram
		    print OUTHISTOGRAM "\t$qtag_best_hist[$index]";
		}
		print OUTHISTOGRAM "\n";
		print OUTHISTOGRAM "$qtag:$stag:second";
		foreach my $index (0..100) { # print qtag second histogram
		    print OUTHISTOGRAM "\t$qtag_second_hist[$index]";
		}
		print OUTHISTOGRAM "\n";
		print OUTHISTOGRAM "$stag:$qtag:good_CGN";
		foreach my $index (0..100) { # print stag good ortho histogram
		    print OUTHISTOGRAM "\t$stag_good_ortho_hist[$index]";
		}
		print OUTHISTOGRAM "\n";
		print OUTHISTOGRAM "$stag:$qtag:second_CGN";
		foreach my $index (0..100) { # print stag good ortho histogram
		    print OUTHISTOGRAM "\t$stag_good_second_hist[$index]";
		}
		print OUTHISTOGRAM "\n";
		print OUTHISTOGRAM "$stag:$qtag:best";
		foreach my $index (0..100) { # print stag best histogram
		    print OUTHISTOGRAM "\t$stag_best_hist[$index]";
		}
		print OUTHISTOGRAM "\n";
		print OUTHISTOGRAM "$stag:$qtag:second";
		foreach my $index (0..100) { # print stag second histogram
		    print OUTHISTOGRAM "\t$stag_second_hist[$index]";
		}
		print OUTHISTOGRAM "\n";
	    }

	    my $max_good_cutoff = 1.0;
	    my $max_value = 0;
	    my $cur_value = 0;
	    if ((($good_qids + $good_sids) == 0) || ((($num_qids + $num_sids) / ($good_qids + $good_sids)) >= 10.0)) {
		print STDERR "Warning number of conserved neighborhoods for ($qtag,$stag) ($good_qids,$good_sids) is much less than total ($num_qids,$num_sids)\n";
		$max_good_cutoff = 0.0;
	    } else {
		for (my $index = 100; $index > 0; $index--) {
		    $cur_value += $qtag_good_ortho_hist[$index];
		    $cur_value += $stag_good_ortho_hist[$index];
		    $cur_value -= 3 * $qtag_good_second_hist[$index]; # 3 is somewhat arbitrary
		    $cur_value -= 3 * $stag_good_second_hist[$index];
		    # $cur_value -= $self_hist{$qtag}[$index]; # seems to penalize transposons too much
		    # $cur_value -= $self_hist{$stag}[$index];
		    if ($cur_value > $max_value) {
			$max_value = $cur_value;
			$max_good_cutoff = $index / 100.0;
		    }
		}
	    }

	    my $max_cutoff = 1.0;
	    $max_value = 0;
	    $cur_value = 0;
	    for (my $index = 100; $index > 0; $index--) {
		$cur_value += $qtag_best_hist[$index];
		$cur_value += $stag_best_hist[$index];
		$cur_value -= 3 * $qtag_second_hist[$index]; # 3 is somewhat arbitrary
		$cur_value -= 3 * $stag_second_hist[$index];
		# $cur_value -= $self_hist{$qtag}[$index]; # seems to penalize transposons too much
		# $cur_value -= $self_hist{$stag}[$index];
		if ($cur_value > $max_value) {
		    $max_value = $cur_value;
		    $max_cutoff = $index / 100.0;
		}
	    }

	    if ($max_cutoff < $max_good_cutoff) { #use the larger cutoff
		$max_cutoff = $max_good_cutoff;
	    }
	    if ($max_cutoff > 0.95) { #never use a cutoff of more than 95%
		$max_cutoff = 0.95;
	    }
	    if ($strict_orthos) {
		$paralog_cutoff[$qTagIndex][$sTagIndex] = $paralog_cutoff[$sTagIndex][$qTagIndex] = $max_cutoff;
	    }
	    if ($good_qids > 0) { #do not divide by zero when there are no matches
		$ave_per_id[$qTagIndex][$sTagIndex] = $qcum_per_id / $good_qids;
		$mean_max_BSR[$qTagIndex][$sTagIndex] = $qcum_max_BSR / $good_qids;
	    } else {
		$ave_per_id[$qTagIndex][$sTagIndex] = 0;
		$mean_max_BSR[$qTagIndex][$sTagIndex] = 0;
	    }
	    if ($good_sids > 0) { #do not divide by zero when there are no matches
		$ave_per_id[$sTagIndex][$qTagIndex] = $scum_per_id / $good_sids;
		$mean_max_BSR[$sTagIndex][$qTagIndex] = $scum_max_BSR / $good_sids;
	    } else {
		$ave_per_id[$sTagIndex][$qTagIndex] = 0;
		$mean_max_BSR[$sTagIndex][$qTagIndex] = 0;
	    }
	    $ave_per_id[$qTagIndex][$sTagIndex] = $ave_per_id[$sTagIndex][$qTagIndex] = ($ave_per_id[$qTagIndex][$sTagIndex] + $ave_per_id[$sTagIndex][$qTagIndex]) / 2;
	    $mean_max_BSR[$qTagIndex][$sTagIndex] = $mean_max_BSR[$sTagIndex][$qTagIndex] = ($mean_max_BSR[$qTagIndex][$sTagIndex] + $mean_max_BSR[$sTagIndex][$qTagIndex]) / 2;
	}
    }
    #set paralog cutoffs for proteins in the same genome
    foreach my $qtag (@tag_array)  {  # start looping through genomes by order in tag file
	my $qTagIndex = $TagIndex{$qtag};
	if ($strict_orthos) {
	    my $min_cutoff = 1.0; # not clear if this should be max or min or something else
	    foreach my $stag (@tag_array) { # loop through all other genomes
		my $sTagIndex = $TagIndex{$stag};
		if ($qtag eq $stag) { # no value stored for same genome cutoffs yet
		    next;
		}
		if ($paralog_cutoff[$qTagIndex][$sTagIndex] < $min_cutoff) {
		    $min_cutoff = $paralog_cutoff[$qTagIndex][$sTagIndex];
		}
	    }
	    $paralog_cutoff[$qTagIndex][$qTagIndex] = $min_cutoff;
	} else {
	    $paralog_cutoff[$qTagIndex][$qTagIndex] = 0.95;
	}
	$ave_per_id[$qTagIndex][$qTagIndex] = 100;
	$mean_max_BSR[$qTagIndex][$qTagIndex] = 1;
    }

    close (OUTHISTOGRAM) if ($histogramfile);

    print STDERR "BSR * 100 Paralog pairwise cutoffs\n";
    open (my $pcutofffile, ">$outprefix" . "pairwise_cutoffs_matrix.txt");
    &print_pairwise_matrix (\@paralog_cutoff, $pcutofffile, 100, 0, 0, 0);
    close ($pcutofffile);

    print STDERR "Mean percent identity for high quality pairwise Orthologs\n";
    open (my $pidentfile, ">$outprefix" . "pairwise_identity_matrix.txt");
    &print_pairwise_matrix (\@ave_per_id, $pidentfile, 1, 0, 0, 0);
    close ($pidentfile);

    print STDERR "Mean maximuim BSR (norm_score) for high quality pairwise Orthologs\n";
    open (my $pBSRfile, ">$outprefix" . "pairwise_BSR_matrix.txt");
    &print_pairwise_matrix (\@mean_max_BSR, $pBSRfile, 100, 0, 0, 0);
    close ($pBSRfile);

    print STDERR "Mean maximuim BSR (norm_score) for high quality pairwise Orthologs Phylip format(distance)\n";
    open (my $pBSRdistpfile, ">$outprefix" . "pairwise_BSR_distance_matrix_phylip.txt");
    &print_pairwise_matrix (\@mean_max_BSR, $pBSRdistpfile, -100, 100, 0, 1);
    close ($pBSRdistpfile);

    print STDERR "Mean maximuim BSR (norm_score) for high quality pairwise Orthologs (distance)\n";
    open (my $pBSRdistfile, ">$outprefix" . "pairwise_BSR_distance_matrix.txt");
    &print_pairwise_matrix (\@mean_max_BSR, $pBSRdistfile, -100, 100, 0, 0);
    close ($pBSRdistfile);

}

sub calc_synbibest  { # subroutine to find bidrectional best synteny scores and populate Qbytaghash and relationship_hash with them

    my $qid = ""; # query feat_name
    my $sid = ""; # subject feat_name
    my $qtag = ""; # query genome tag
    my $stag = ""; # subject tag

    #mark best synteny score per genome per qid
    foreach $qtag (@tag_array)  {  # start looping through by order in tag file (reference is first)
	foreach $qid ( keys %{ $Tagbyfeatnamehash{$qtag} } )  { # go through featnames of query genome to look for matches in other genomes
	    if (!defined $Qbytaghash{$qid}) {
		next;# need to do this to prevent keys %{ $Qbytaghash{$qid}{$stag} } causing $Qbytaghash{$qid}{$stag} to be defined later
	    }
	    foreach $stag (keys %{ $Qbytaghash{$qid} })  {
		if ($qtag eq $stag) {
		    next; #skip matches within the same genome
		}
		my $best_score = 0;
		my $second_score = 0;
		my $found_best = 0;
		my $best_sid;
		my $second_sid;
		foreach $sid (sort {$Qbytaghash{$qid}{$stag}{$b}->{'synteny'} <=> $Qbytaghash{$qid}{$stag}{$a}->{'synteny'}} keys %{ $Qbytaghash{$qid}{$stag} })  {
		    if (!$found_best) {
			$best_score = $relationship_hash{$qid}{$sid}->{'synteny'};
			$found_best = 1;
			$best_sid = $sid;
			next;
		    }
		    $second_score = $relationship_hash{$qid}{$sid}->{'synteny'};
		    $second_sid = $sid;
		    last; #quit after finding second best synteny score
		}
		if ($found_best) {
		    $relationship_hash{$qid}{$best_sid}->{'synbest'} = 1;
		    if ((0.5 * $best_score) < $second_score) {
			$relationship_hash{$qid}{$second_sid}->{'synbest'} = 1;# allow more than one synbest for realtively good second scores
		    }
		}
	    }
	}
    }

    #mark bidirectionally best synteny scores per genome per qid if it exists
    foreach $qid (keys %relationship_hash)  { # go through all featnames
	foreach $sid (keys %{ $relationship_hash{$qid} } )  {
		if ($relationship_hash{$qid}{$sid}->{'synbest'} && $relationship_hash{$sid}{$qid}->{'synbest'}) {
		    $relationship_hash{$qid}{$sid}->{'synbibest'}  = 1;
		}
	}
    }
			
}

sub calc_cliques  { # subroutine to find bidrectional best cliques and populate Qbytaghash and relationship_hash with them

    foreach my $qid (keys %relationship_hash)  { # go through all featnames
	my %clique_ids = ();
	my %clique_tags = ();
	my $clique_already_tried = 0;
	my $clique_all = 1;
	my $clique_failed = 0;
	my $num_clique_matches;
	my $qtag = $FeatnameLookupTag_hash{$qid};
	my $clique_top = 1;
	foreach my $sid (sort {$relationship_hash{$qid}{$b}->{'BSR'} <=> $relationship_hash{$qid}{$a}->{'BSR'}} keys %{ $relationship_hash{$qid} } )  {
	    print STDERR "Clique $qid $sid ($relationship_hash{$qid}{$sid}->{'BSR'}) :$relationship_hash{$qid}{$sid}->{'bibest'}\n" if ($DEBUG);
	    if ($qid eq $sid) {
		next; #skip matches to same protein
	    }
	    my $stag = $FeatnameLookupTag_hash{$sid};
	    if ($relationship_hash{$qid}{$sid}->{'bibest'}) {
		if (defined $clique_tags{$stag}){
		    print STDERR "WARNING!!! more than one bibest for qid $qid $qtag: stag $stag: $sid\n";
		}
		if (defined $relationship_hash{$qid}{$sid}->{'clique_all'}) {
		    $clique_already_tried = 1;
		    print STDERR "already tried :$relationship_hash{$qid}{$sid}->{'clique_all'}\n" if ($DEBUG);
		    last;
		}
		$clique_ids{$sid} = 1;
		$clique_tags{$stag} = 1;
	    } else {
		$clique_all = 0; #all matches are not part of the clique
		last; #end of possible clique
	    }
	}
	if ($clique_already_tried) {
	    next; #tried this clique already from some other starting point
	}
	$num_clique_matches = keys %clique_ids;
	if ($num_clique_matches < 2) {
	    $clique_failed = 1;
	    $clique_ids{$qid} = 1; #add the query into the clique
	} else {
	    my @clique_ids_to_check = keys %clique_ids;
	    $clique_ids{$qid} = 1; #add the query into the clique
	    foreach my $cqid (@clique_ids_to_check) {
		my $num_cqid_clique_matches = 0;
		my $clique_cqid_failed = 0;
		my @sort_array = sort {$relationship_hash{$cqid}{$b}->{'BSR'} <=> $relationship_hash{$cqid}{$a}->{'BSR'}} keys %{ $relationship_hash{$cqid} };
		my $num_cqid_matches = @sort_array;
		foreach my $csid (@sort_array)  {
		    if ($cqid eq $csid) {
			$num_cqid_matches--;
			next; #skip matches to same protein
		    }
		    if ($relationship_hash{$cqid}{$csid}->{'bibest'}) {
			if (!defined $clique_ids{$csid}){
			    $clique_cqid_failed = 1; #should not have an additional bibest not in clique
			    last;
			}
			$num_cqid_clique_matches++;
			if (defined $relationship_hash{$cqid}{$csid}->{'clique_all'}) {
			    print STDERR "WARNING!!! this clique already marked as tried for qid $qid $qtag: cqid $cqid: csid $csid\n" if ($DEBUG);
			    last;
			}
		    } else {
			last; #end of possible clique
		    }
		}
		if (($clique_cqid_failed) || ($num_cqid_clique_matches != $num_clique_matches)) {
		    $clique_failed = 1;
		    last;
		}
		if ($num_cqid_matches != $num_clique_matches) {
		    $clique_all = 0;
		}
	    }
	}
	if ($clique_failed) {
	    $clique_all = 0;
	    $clique_top = 0;
	} else {
	    $clique_top = $num_clique_matches + 1; # include the query in the clique count
	    if ($clique_all) {
		$clique_all = $num_clique_matches + 1; # include the query in the clique count
	    }
	}
	if ($DEBUG){
	    my @tmp_clique_ids = (keys %clique_ids);
	    print STDERR "$clique_all : $clique_top @tmp_clique_ids\n";
	}
	foreach my $cqid (keys %clique_ids) {
	    foreach my $csid (keys %clique_ids) {
		if ($cqid eq $csid) {
		    next; #skip matches to same protein
		}
		if (defined $relationship_hash{$cqid}{$csid}) {
		    $relationship_hash{$cqid}{$csid}->{'clique_all'} = $clique_all;
		    $relationship_hash{$cqid}{$csid}->{'clique_top'} = $clique_top;
		}
	    }
	}
    }
			
    #mark any unmarked clique_top and clique_all as 0
    foreach my $qid (keys %relationship_hash)  { # go through all featnames
	foreach my $sid (keys %{ $relationship_hash{$qid} } )  {
	    if (!defined $relationship_hash{$qid}{$sid}->{'clique_top'}) {
		$relationship_hash{$qid}{$sid}->{'clique_all'} = 0;
		$relationship_hash{$qid}{$sid}->{'clique_top'} = 0;
	    }
	}
    }
			
}

sub calc_synteny  { # subroutine to calculate synteny scores and populate Qbytaghash and relationship_hash with them

    my ($window) = @_;

   #call synteny score for each relationship
    foreach my $qid (keys %relationship_hash)  { # go through all featnames
	foreach my $sid (keys %{ $relationship_hash{$qid} } )  {
	    if ($qid eq $sid) {
		next; #skip matches to same protein
	    }
	    my $qtag = $FeatnameLookupTag_hash{$qid};
	    if (!defined $qtag) {
		print STDERR "$qid has no defined tag in FeatnameLookupTag_hash\n";
	    }
	    my $stag = $FeatnameLookupTag_hash{$sid};
	    if (!defined $stag) {
		print STDERR "$sid has no defined tag in FeatnameLookupTag_hash\n";
	    }
	    if ($qtag eq $stag) {
		next; #skip matches to same genome
	    }
	    print STDERR "$qid $sid $relationship_hash{$qid}{$sid}->{'id'} $relationship_hash{$qid}{$sid}->{'eval'} $relationship_hash{$qid}{$sid}->{'score'} $relationship_hash{$qid}{$sid}->{'best'} $relationship_hash{$qid}{$sid}->{'bibest'} $relationship_hash{$qid}{$sid}->{'clique_top'} $relationship_hash{$qid}{$sid}->{'clique_all'}\n" if($DEBUG);
	    ($relationship_hash{$qid}{$sid}->{'synteny'}, $relationship_hash{$qid}{$sid}->{'CGN_bibest'}) = &synteny($qtag, $qid, $stag, $sid, $window);
	}
    }
			
}

sub merge_clusters {#&merge_clusters(qid, sid) is used to test if two clusters can be merged and if so to merge them
    my ($qid, $sid) = @_;
    my $merged = [];
    my $failed_merge = 0;
    
    #assumes the cluster arrays are sorted by genome tag index
    my $index = 0; #need to keep track of where we are in the subject array
    foreach my $cur_qid_tag ( @{ $clusters{$qid} } ) {
	while (($index <= $#{ $clusters{$sid} }) && (($clusters{$sid}->[$index])->{'tag'} < $cur_qid_tag->{'tag'})) {
	    push(@{ $merged }, $clusters{$sid}->[$index]);
	    print STDERR Dumper($merged) if ($DEBUG);
	    $index++;
	}
	if ($index <= $#{ $clusters{$sid} }) {
	    if (($clusters{$sid}->[$index])->{'tag'} == $cur_qid_tag->{'tag'}) {
		if (($clusters{$sid}->[$index])->{'id'} eq $cur_qid_tag->{'id'}) {
		    die ("ERROR: two different clusters should never contain the same protein\n");
		} else {
		    $failed_merge = 1; #cannot have two different proteins from the same genome
		    last;
		}
	    }
	}
	push(@{ $merged }, $cur_qid_tag);
	print STDERR Dumper($merged) if ($DEBUG);
    }
    if (!$failed_merge) {
	while ($index <= $#{ $clusters{$sid} }) {
	    push(@{ $merged }, $clusters{$sid}->[$index]);
	    print STDERR Dumper($merged) if ($DEBUG);
	    $index++;
	}
	foreach my $cur_id_tag ( @{ $clusters{$qid} }, @{ $clusters{$sid} } ) {
	    $clusters{$cur_id_tag->{'id'}} = $merged;
	}
    } else {
	print STDERR "failed merge\n" if ($DEBUG);
    }
    
}

sub calc_clusters { # greedily compute clusters by starting with largest relationship score to merge existing clusters (starting with every protein as a singleton) and constrainging a cluster to only contain one protein from each genome

    # start by creating a sorted array of relationships from the relationship_hash
    my @relationship_array = ();
    foreach my $qid (keys %relationship_hash)  { # go through all featnames
	foreach my $sid (keys %{ $relationship_hash{$qid} } )  { # go through all relationships
	    if ($qid eq $sid) {
		next; #skip matches to same protein
	    }
	    my $qtag = $FeatnameLookupTag_hash{$qid};
	    my $qTagIndex = $TagIndex{$qtag};
	    my $stag = $FeatnameLookupTag_hash{$sid};
	    my $sTagIndex = $TagIndex{$stag};
	    if ($qtag eq $stag) {
		next; #skip matches to same genome
	    }
	    if (!$relationship_hash{$qid}{$sid}->{'synbibest'}) {
		next; # ignore matches which are not bidirectionally best synteny matches
	    }
	    if (($strict_orthos == 1) && !($relationship_hash{$qid}{$sid}->{'anchor'} || $relationship_hash{$qid}{$sid}->{'extend'} || ($relationship_hash{$qid}{$sid}->{'CGN_bibest'} > $CGN_window_size) || ((($relationship_hash{$qid}{$sid}->{'norm_score'} >= $paralog_cutoff[$qTagIndex][$sTagIndex]) || (defined $relationship_hash{$sid}{$qid} && ($relationship_hash{$sid}{$qid}->{'norm_score'} >= $paralog_cutoff[$sTagIndex][$qTagIndex]))) && ($relationship_hash{$qid}{$sid}->{'bibest'} || ($relationship_hash{$qid}{$sid}->{'CGN_bibest'} > 1))))) {
		next; # ignore matches which do not meet the strict ortholog criteria if that option was requested
		# strict criteria are having at least 1/2 of the maximal number of possible bibest CGN matches or
		# having a BSR above the paralog cutoff and either being a bibest match or having 1 or more CGN bibest matches
	    }
	    if (($strict_orthos == 2) && !($relationship_hash{$qid}{$sid}->{'anchor'} || $relationship_hash{$qid}{$sid}->{'extend'} || ($relationship_hash{$qid}{$sid}->{'CGN_bibest'} > $CGN_window_size))) {
		next; # ignore matches which do not meet the strict ortholog criteria if that option was requested
		# strict criteria are having at least 1/2 of the maximal number of possible bibest CGN matches or
	    }
	    my $score = $relationship_hash{$qid}{$sid}->{'synteny'};
	    if (defined $relationship_hash{$sid}{$qid}) {
		if ($qid gt $sid) {
		    next; # only enter qid,sid pairs once not symmetrically since they indicate the same cluster join
		}
		if ($relationship_hash{$sid}{$qid}->{'synteny'} > $score) {# use the larger synteny score
		    $score = $relationship_hash{$sid}{$qid}->{'synteny'};
		}
	    }
	    my $hash_ref = {};
	    $hash_ref->{'query'} = $qid;
	    $hash_ref->{'subject'} = $sid;
	    $hash_ref->{'score'} = $score;
	    push (@relationship_array, $hash_ref);
	}
    }
    @relationship_array = sort {$b->{'score'} <=> $a->{'score'}} ( @relationship_array );
    print STDERR Dumper(\@relationship_array) if ($DEBUG);

    foreach my $qid (keys %relationship_hash)  { # go through all featnames and initialize clusters as singletons
	if (!defined $FeatnameLookupTag_hash{$qid}) {
	    die("ERROR!!! $qid in relationship_hash but not in FeatnameLookupTag_hash\n");
	}
	$clusters{$qid} = [{ 'id' => $qid, 'tag' => $TagIndex{$FeatnameLookupTag_hash{$qid}} }];
    }
    print STDERR Dumper(\%clusters) if ($DEBUG);
    foreach my $match ( @relationship_array ) {
	if ($match->{'score'} < $min_synteny_threshold) {
	    last; # stop using matches when the scores get too low
	}
	if ($clusters{$match->{'query'}} == $clusters{$match->{'subject'}}) {
	    next; #already in the same cluster
	}
	print STDERR "merge ($match->{'query'}, $match->{'subject'}): $match->{'score'}\n" if ($DEBUG);
	&merge_clusters($match->{'query'}, $match->{'subject'});
    }
    print STDERR "Done clustering\n" if ($DEBUG);
    print STDERR Dumper(\%clusters) if ($DEBUG);
			
}

sub calc_cluster_numbers { # order the cluster numbers by the reference genome(s)

############ Assign cluster numbers to clusters in the same order they will be output in the match tables and compute cluster_size ###############
    my $cluster_num = 0;
    
    print STDERR "Calculating cluster numbers ...\n";
    foreach my $genome_tag (@tag_array)  {  # use order in tag file (reference is first) so clusters are numbered in order we want them output
	foreach my $asmbl_id (sort keys %{ $genome_hash{$genome_tag} }) { #output by ascending asmbl_id
	    for my $index (0 .. $#{ $genome_hash{$genome_tag}{$asmbl_id} }) { #output by ascending 5'(3') protein coordinate
		if (defined $genome_hash_context{$genome_tag}{$asmbl_id}->{$index}) {
		    next; #skip features/genes/proteins that are just there for context
		}
		my $query_featname = $genome_hash{$genome_tag}{$asmbl_id}->[$index];
		$cluster_num++;
		$cluster_size[$cluster_num] = 0;
		foreach my $cluster_member ( @{ $clusters{$query_featname} } ) {
		    if (defined $cluster_number{$cluster_member->{'id'}}) { # we have seen this cluster before so skip and reset counter
			$cluster_num--;
			last;
		    } else {
			$cluster_number{$cluster_member->{'id'}} = $cluster_num;
			$cluster_for_number[$cluster_num] = $clusters{$query_featname};
			$cluster_size[$cluster_num]++;
		    }
		}
	    }
	}
    }
    $#cluster_size = $cluster_num; # need to make sure we only keep real values since we sometimes deprecate a prospective cluster (see $cluster_num-- above)
    return;
}

sub calc_cluster_weights { # record strong and weak matches between clusters and very strong matches within clusters

    my %paralog_weight = (); # Key1 = cluster_number_1, Key2 = cluster_number_2, Value = number of strong matches between clusters
    my %cluster_weight = (); # Key = cluster number, Value = number of very strong matches within a cluster
    my @cluster_length = (); # Array of hashs, Index = cluster number, Hash stores min, max, mean, std dev

    foreach my $cluster_num (1 .. $#cluster_size)  {
	my $max_BSR_sum = -1;
	my $centroid = undef;
	my $mymin = 1000000;
	my $mymax = 0;
	my $mymean = 0;
	my $mymeansquared = 0;
	$cluster_weight{$cluster_num} = 0; # set cluster_weight to 0 here for convenience
	print STDERR "cluster $cluster_num\n" if ($DEBUG);
	foreach my $cluster_member ( @{ $cluster_for_number[$cluster_num] } ) {
	    my $protein_length = $feat_hash{$cluster_member->{'id'}}->{'length'};
	    if ($mymin > $protein_length) {
		$mymin = $protein_length;
	    }
	    if ($mymax < $protein_length) {
		$mymax = $protein_length;
	    }
	    $mymean += $protein_length;
	    $mymeansquared += $protein_length * $protein_length;
	}
	$cluster_length[$cluster_num] = {};
	$cluster_length[$cluster_num]{'min'} = $mymin;
	$cluster_length[$cluster_num]{'max'} = $mymax;
	$mymean /= $cluster_size[$cluster_num];
	$cluster_length[$cluster_num]{'mean'} = $mymean;
	my $normalized_cluster_size = $cluster_size[$cluster_num];
	if ($normalized_cluster_size > 1) {
	    $normalized_cluster_size--;
	}
	$cluster_length[$cluster_num]{'stddev'} = sqrt(($mymeansquared - ($mymean * $mymean * $cluster_size[$cluster_num])) / $normalized_cluster_size);
    }

    if ($fragmentfile)  {
	print STDERR "Finding possible protein fragments and fusions...\n";
	open (OUTFRAGFUS, ">$outprefix" . "fragments_fusions.txt");
    }
    if ($strict_orthos)  {
	print STDERR "Finding below cutoff cluster members...\n";
	open (OUTSTRICT, ">$outprefix" . "below_cutoff_clusters.txt");
    }

    foreach my $qid (keys %relationship_hash)  { # go through all featnames
	my $matches_in_cluster = 0;
	my $matches_out_of_cluster = 0;
	my $matches_above_cutoff = 0;
	my $full_length = 0;
	my $short_paralog = 0;
	my $qtag = $FeatnameLookupTag_hash{$qid};
	my $qTagIndex = $TagIndex{$qtag};
	my $cluster_number_qid = $cluster_number{$qid};
	foreach my $sid (keys %{ $relationship_hash{$qid} } )  { # go through all relationships
	    if ($qid eq $sid) {
		next; #skip matches to same protein
	    }
	    if (!defined $cluster_number{$qid}) {
		die ("cluster number not defined for $qid\n");
	    }
	    if (!defined $cluster_number{$sid}) {
		die ("cluster number not defined for $sid\n");
	    }
	    my $stag = $FeatnameLookupTag_hash{$sid};
	    my $sTagIndex = $TagIndex{$stag};
	    my $cluster_number_sid = $cluster_number{$sid};
	    if ($cluster_number_qid == $cluster_number_sid) {
		$matches_in_cluster++;
		if ($relationship_hash{$qid}{$sid}->{'full'}) {
		    $full_length++;
		}
		if ($strict_orthos && (($relationship_hash{$qid}{$sid}->{'norm_score'} >= $paralog_cutoff[$qTagIndex][$sTagIndex]) || (defined $relationship_hash{$sid}{$qid} && ($relationship_hash{$sid}{$qid}->{'norm_score'} >= $paralog_cutoff[$sTagIndex][$qTagIndex])))) {
		    $matches_above_cutoff++;
		}
		if ($qid gt $sid) {
		    next; # only enter qid,sid pairs once not symmetrically since they indicate the same cluster weight
		}
		if (($relationship_hash{$qid}{$sid}->{'synbibest'}) && (!$strict_orthos || $relationship_hash{$qid}{$sid}->{'anchor'} || $relationship_hash{$qid}{$sid}->{'extend'} || ($relationship_hash{$qid}{$sid}->{'CGN_bibest'} > $CGN_window_size) || ((($relationship_hash{$qid}{$sid}->{'norm_score'} >= $paralog_cutoff[$qTagIndex][$sTagIndex]) || (defined $relationship_hash{$sid}{$qid} && ($relationship_hash{$sid}{$qid}->{'norm_score'} >= $paralog_cutoff[$sTagIndex][$qTagIndex]))) && ($relationship_hash{$qid}{$sid}->{'bibest'} || ($relationship_hash{$qid}{$sid}->{'CGN_bibest'} > 1))))) {
		    if (!defined $cluster_weight{$cluster_number_qid}) {
			die ("all cluster_weights should have been set to zero when clusters were assigned numbers\n");
		    } else {
			$cluster_weight{$cluster_number_qid}++;
		    }
		}
		next;
	    }
	    if (!$strict_orthos || $relationship_hash{$qid}{$sid}->{'anchor'} || $relationship_hash{$qid}{$sid}->{'extend'} || ($relationship_hash{$qid}{$sid}->{'CGN_bibest'} > $CGN_window_size) || (($relationship_hash{$qid}{$sid}->{'norm_score'} >= $paralog_cutoff[$qTagIndex][$sTagIndex]) || (defined $relationship_hash{$sid}{$qid} && ($relationship_hash{$sid}{$qid}->{'norm_score'} >= $paralog_cutoff[$sTagIndex][$qTagIndex])))) {
		$matches_out_of_cluster++;
		if (!$relationship_hash{$qid}{$sid}->{'full'}) {
		    $short_paralog++;
		}
		if ($qid gt $sid) {
		    next; # only enter qid,sid pairs once not symmetrically since they indicate the same cluster weight
		}
		my $first_cluster_number = $cluster_number_qid;
		my $second_cluster_number = $cluster_number_sid;
		if ($cluster_number_qid > $cluster_number_sid) { # store in cluster weight hashes only for $cluster_number_qid < $cluster_number_sid
		    $first_cluster_number = $cluster_number_sid;
		    $second_cluster_number = $cluster_number_qid;
		}
		if (!defined $paralog_weight{$first_cluster_number}{$second_cluster_number}) {
		    $paralog_weight{$first_cluster_number}{$second_cluster_number} = 1;
		} else {
		    $paralog_weight{$first_cluster_number}{$second_cluster_number}++;
		}
	    }
	}
	print STDERR "$qid $cluster_number_qid $feat_hash{$qid}->{'retained'} $feat_hash{$qid}->{'length'} $cluster_length[$cluster_number_qid]->{'mean'} $full_length $matches_in_cluster $short_paralog $matches_out_of_cluster\n" if ($DEBUG);
	if ($fragmentfile && !$feat_hash{$qid}->{'retained'}) {#we have already output frameshift fragments so don't do it again
	    if (($feat_hash{$qid}->{'length'} < (0.9 * $cluster_length[$cluster_number_qid]->{'mean'})) && ((2 * $full_length) < $matches_in_cluster)) {
		print OUTFRAGFUS "$qid\tFragment in cluster\n";
	    } elsif (($feat_hash{$qid}->{'length'} > (1.1 * $cluster_length[$cluster_number_qid]->{'mean'})) && ($short_paralog > $full_length) && ($short_paralog == $matches_out_of_cluster) && (($short_paralog + $full_length + 1) >= $matches_in_cluster)) {
		print OUTFRAGFUS "$qid\tFusion\n";
	    } elsif (($matches_out_of_cluster > $matches_in_cluster) && ($short_paralog == $matches_out_of_cluster)) {
		print OUTFRAGFUS "$qid\tFragment out of cluster\n";
	    }
	}
	if ($strict_orthos && ((2 * $matches_above_cutoff) < $matches_in_cluster)) {
	    print OUTSTRICT "$cluster_number_qid\t$qid\n";
	}
    }
    print STDERR "Done cluster weights\n" if ($DEBUG);
    if ($fragmentfile)  {
	close (OUTFRAGFUS);
    }
    if ($strict_orthos)  {
	close (OUTSTRICT);
    }
			
    if ($writeparalogs)  {
	print STDERR "Finding paralogs...\n";

	open (OUTPARAW, ">$outprefix" . "paralog_weights.txt");
	foreach my $cluster_id1 (sort { $a <=> $b } keys %paralog_weight)  {  # loop through all clusters
	    foreach my $cluster_id2 ( sort { $a <=> $b } keys %{ $paralog_weight{$cluster_id1} } )  { # by all clusters with weights
		print OUTPARAW "$cluster_id1\t$cluster_id2\t$paralog_weight{$cluster_id1}{$cluster_id2}\n";
	    }
	}
	close (OUTPARAW);

	open (OUTPARALOGS, ">$outprefix" . "paralogs.txt");
	my %used_in_paralog = ();
	# make paralog_weight symmetrical for paralog clustering
	foreach my $cluster_id1 (sort { $a <=> $b } keys %paralog_weight)  {  # loop through all clusters
	    foreach my $cluster_id2 ( sort { $a <=> $b } keys %{ $paralog_weight{$cluster_id1} } )  { # by all clusters with weights
		$paralog_weight{$cluster_id2}{$cluster_id1} = $paralog_weight{$cluster_id1}{$cluster_id2};
	    }
	}
	foreach my $cluster_id1 (sort { $a <=> $b } keys %paralog_weight)  {  # loop through all clusters
	    if (defined $used_in_paralog{$cluster_id1}) {
		next; # already included this in a paralog cluster
	    }
	    my $first = 1;
	    my @paralog_array = ();
	    push (@paralog_array, $cluster_id1);
	    my $next_cluster_id;
	    while (defined ($next_cluster_id = shift(@paralog_array))) {
		foreach my $cluster_id2 ( sort { $a <=> $b } keys %{ $paralog_weight{$next_cluster_id} } )  { # by all clusters with weights
		    if (defined $used_in_paralog{$cluster_id2}) {
			next; # already included this in a paralog cluster
		    }
		    if ((!defined $relationship_hash{$centroids[$cluster_id1]}{$centroids[$cluster_id2]}) || (!$relationship_hash{$centroids[$cluster_id1]}{$centroids[$cluster_id2]}->{'full'})) {
			next; # only call paralogs which are full length matches not domain matches
		    }
		    $used_in_paralog{$cluster_id2} = 1;
		    if ($first) {
			print OUTPARALOGS "$cluster_id1";
			$used_in_paralog{$cluster_id1} = 1;
			$first = 0;
		    }
		    print OUTPARALOGS "\t$cluster_id2";
		    push (@paralog_array, $cluster_id2);
		}
	    }
	    if (!$first) {
		print OUTPARALOGS "\n";
	    }
	}
	close (OUTPARALOGS);
	open (OUTPARALOGS, ">$outprefix" . "partial_paralogs.txt");
	#output non full length paralogs which did not get grouped with other clusters
	foreach my $cluster_id1 (sort { $a <=> $b } keys %paralog_weight)  {  # loop through all clusters
	    if (defined $used_in_paralog{$cluster_id1}) {
		next; # already included this in a paralog cluster
	    }
	    print OUTPARALOGS "$cluster_id1";
	    foreach my $cluster_id2 ( sort { $a <=> $b } keys %{ $paralog_weight{$cluster_id1} } )  { # by all clusters with weights
		print OUTPARALOGS "\t$cluster_id2";
	    }
	    print OUTPARALOGS "\n";
	}
	close (OUTPARALOGS);
    }

    if ($writeclusterweights)  {
	print STDERR "Finding cluster weights...\n";
	open (OUTCWEIGHT, ">$outprefix" . "cluster_weights.txt");
	foreach my $cluster_num (sort { $a <=> $b } keys %cluster_weight)  {  # loop through all clusters
	    my $max_weight = ($cluster_size[$cluster_num] * ($cluster_size[$cluster_num] - 1)) / 2;
	    my $per_weight = 100;
	    if ($max_weight > 0) {
		$per_weight = 100 * ($cluster_weight{$cluster_num} / $max_weight);
	    }
	    printf OUTCWEIGHT "%d\t%d\t%5.2f\t%d\t%d\t%7.2f\t%7.2f\n", $cluster_num, $cluster_size[$cluster_num], $per_weight, $cluster_length[$cluster_num]->{'min'}, $cluster_length[$cluster_num]->{'max'}, $cluster_length[$cluster_num]->{'mean'}, $cluster_length[$cluster_num]->{'stddev'};
	}
	close (OUTCWEIGHT);
    }

}

sub calc_anchors { #calculate anchors which are matches with $ANCHOR_WINDOW bibest matches on both sides of the match

    foreach my $qid (keys %relationship_hash)  { # go through all featnames
	foreach my $sid (keys %{ $relationship_hash{$qid} } )  { # go through all relationships
	    if ($qid eq $sid) {
		next; #skip matches to same protein
	    }
	    my $qtag = $FeatnameLookupTag_hash{$qid};
	    my $stag = $FeatnameLookupTag_hash{$sid};
	    if ($qtag eq $stag) {
		next; #skip matches to same genome
	    }
	    if ($relationship_hash{$qid}{$sid}->{'CGN_bibest'} >= $ANCHOR_cutoff) {
		$relationship_hash{$qid}{$sid}->{'anchor'} = 1;
	    }
	}
    }
}

sub extend_anchors { #extend anchors by marching out from anchors

    foreach my $qid (keys %relationship_hash)  { # go through all featnames
	foreach my $sid (keys %{ $relationship_hash{$qid} } )  { # go through all relationships
	    if ($qid eq $sid) {
		next; #skip matches to same protein
	    }
	    my $qtag = $FeatnameLookupTag_hash{$qid};
	    my $stag = $FeatnameLookupTag_hash{$sid};
	    if ($qtag eq $stag) {
		next; #skip matches to same genome
	    }
	    if (!$relationship_hash{$qid}{$sid}->{'anchor'}) {
		print STDERR "$qid $sid not an anchor\n" if ($DEBUG);
		next; #not an anchor so skip
	    }
	    my $QueryArrayIndex = $TagByPointer{$qid}; # points to the position in the genome_hash array
	    my $SubjectArrayIndex = $TagByPointer{$sid}; # points to the position in the genome_hash array
	    my $query_orient = $feat_hash{$qid}->{'orient'};
	    my $subject_orient = $feat_hash{$sid}->{'orient'};
	    my $orient = $query_orient * $subject_orient;
	    my $inc_dec = "";
	    my $stop_j = "";
	    my $i = "";
	    my $j = "";
	    my $query = "";
	    my $subject = "";
	    my $Qasmbl_id = $AssemblyLookup_hash{$qid};
	    my $Sasmbl_id = $AssemblyLookup_hash{$sid};    
	    my $query_array_last = $#{ $genome_hash{$qtag}{$Qasmbl_id} };
	    my $subject_array_last = $#{ $genome_hash{$stag}{$Sasmbl_id} };
	    my $skip_tries = 0;

	    print STDERR "query $qid ($QueryArrayIndex : $query_array_last)     subject $sid ($SubjectArrayIndex : $subject_array_last)\n" if ($DEBUG);

	    $j = $SubjectArrayIndex;
	    if ($orient > 0) {
		$inc_dec = 1; #anchor proteins are on the same strand
		$stop_j = $subject_array_last;
	    } else {
		$inc_dec = -1; #anchor proteins are on opposite strands
		$stop_j = 0;
	    }
	    for ($i = $QueryArrayIndex + 1; $i <= $query_array_last; $i++) { # iterate upward first
		$query = ${ $genome_hash{$qtag}{$Qasmbl_id} }[$i];
		print STDERR "       qpos = $i ($query)\n" if ($DEBUG);
		if ($j == $stop_j) {
		    print STDERR "         end of assembly\n" if ($DEBUG);
		    last; #ran off the end of the assembly
		}
		$j += $inc_dec;
		$subject = ${ $genome_hash{$stag}{$Sasmbl_id} }[$j];
		if (!defined $relationship_hash{$query}{$subject}) {
		    print STDERR "         $subject $j no match\n" if ($DEBUG);
		    if ($skip_tries < 2) {
			$skip_tries++;
			print STDERR "         skipping $query $i\n" if ($DEBUG);
			$j -= $inc_dec; #relationship_hash is symmetric so skipping should occur on the other genome as well but only once - cannot mix
			next;
		    } else {
			last; #exceeded skipping over protein tries - no match so cannot extend
		    }
		}
		if ($relationship_hash{$query}{$subject}->{'anchor'}) {
		    print STDERR "         $subject $j already an anchor\n" if ($DEBUG);
		    last; #already an anchor so can extend from it instead
		}
		print STDERR "         spos = $j ($subject) [$relationship_hash{$query}{$subject}->{'id'}] <=> [$relationship_hash{$subject}{$query}->{'id'}] $relationship_hash{$query}{$subject}->{'best'} $relationship_hash{$query}{$subject}->{'bibest'} $relationship_hash{$query}{$subject}->{'clique_top'} $relationship_hash{$query}{$subject}->{'clique_all'} (orient: $orient $feat_hash{$query}->{'orient'} $feat_hash{$subject}->{'orient'})\n" if ($DEBUG);
		$relationship_hash{$query}{$subject}->{'extend'} = 1;
		$skip_tries = 0; #skip tries are not cummulative can repeatedly skip two proteins after a match
	    }

	    print STDERR "     going down now\n" if ($DEBUG);

	    $j = $SubjectArrayIndex;
	    if ($orient > 0) {
		$inc_dec = -1; #anchor proteins are on the same strand
		$stop_j = 0;
	    } else {
		$inc_dec = 1; #anchor proteins are on opposite strands
		$stop_j = $subject_array_last;
	    }
	    $skip_tries = 0; #reset
	    for ($i = $QueryArrayIndex - 1; $i >= 0; $i--) { # iterate downward now
		$query = ${ $genome_hash{$qtag}{$Qasmbl_id} }[$i];
		print STDERR "       qpos = $i ($query)\n" if ($DEBUG);
		if ($j == $stop_j) {
		    print STDERR "         end of assembly\n" if ($DEBUG);
		    last; #ran off the end of the assembly
		}
		$j += $inc_dec;
		$subject = ${ $genome_hash{$stag}{$Sasmbl_id} }[$j];
		if (!defined $relationship_hash{$query}{$subject}) {
		    print STDERR "         $subject $j no match\n" if ($DEBUG);
		    if ($skip_tries < 2) {
			$skip_tries++;
			print STDERR "         skipping $query $i\n" if ($DEBUG);
			$j -= $inc_dec; #relationship_hash is symmetric so skipping should occur on the other genome as well but only once - cannot mix
			next;
		    } else {
			last; #exceeded skipping over protein tries - no match so cannot extend
		    }
		}
		if ($relationship_hash{$query}{$subject}->{'anchor'}) {
		    print STDERR "         $subject $j already an anchor\n" if ($DEBUG);
		    last; #already an anchor so can extend from it instead
		}
		print STDERR "         spos = $j ($subject) [$relationship_hash{$query}{$subject}->{'id'}] <=> [$relationship_hash{$subject}{$query}->{'id'}] $relationship_hash{$query}{$subject}->{'best'} $relationship_hash{$query}{$subject}->{'bibest'} $relationship_hash{$query}{$subject}->{'clique_top'} $relationship_hash{$query}{$subject}->{'clique_all'} (orient: $orient $feat_hash{$query}->{'orient'} $feat_hash{$subject}->{'orient'})\n" if ($DEBUG);
		$relationship_hash{$query}{$subject}->{'extend'} = 1;
		$skip_tries = 0; #skip tries are not cummulative can repeatedly skip two proteins after a match
	    }
	}
    }
}

sub print_adjacency_matrix { # print the entries of the adjacecny matrix for one node, return the max weight node and strand

    my ($file, $node, $adjacency_matrix_ref) = @_;
    my $first = 1;
    my $max = undef;

    if ($DEBUG) {
	print STDERR "PAMin: $node\n";
    }

    foreach my $next_node (sort {${ $adjacency_matrix_ref }{$node}{$b} <=> ${ $adjacency_matrix_ref }{$node}{$a}} keys %{ ${ $adjacency_matrix_ref }{$node} } ) {
	if ($first) {
	    $first = 0;
	    $max = $next_node;
	} else {
	    print $file "\t";
	}
	print $file "($node,$next_node,${ $adjacency_matrix_ref }{$node}{$next_node})";
    }
    if (defined $max) {
	my ($cluster_id, $end) = split ('_', $max, 2);
	my $strand;
	if ($end eq "5") {
	    $strand = "+";
	} elsif ($end eq "3") {
	    $strand = "-";
	} else {
	    die ("ERROR: $max is an improperly formatted node in the adjacency matrix!\n");
	}
	if ($DEBUG) {
	    print STDERR "PAMout:$max $cluster_id $strand\n";
	}

	return ($cluster_id, $strand);
    } else {
	return (undef, undef);
    }

}

sub calc_adjacency_weights { # construct adjacency matrix for given level of cluster representation (genomes in cluster)

    my ($level) = @_;

    my $size_cutoff = ($level * $genome_number) / 100;
    my %adjacency_matrix = (); # Key1 = cluster_number_1+end(_5 or _3), Key2 = cluster_number_2+end(_5 or_3), Value = number of times clusters are adjacent at the given level
    my %cluster_visited = (); #Key = cluster_number, Value = doesn't matter just being defined marks as visited
    my %hash_breaks = (); #Key = cluster_number_2+end(_5 or_3), Value = number of different rearrangements at this cluster end, only computed for level = 100 where it makes sense
    my @vector_breaks = (); #at each node different rearrangements can occur or an adjacency can be missing due to a contig break, use this vector to track this difference
    my @pairwise_adj_breaks = (); #matrix of number of adjacency breaks between pairs of genomes, keep symmetric, output as a distance matirx
    
    for ( my $index1 = 0 ; $index1 < $genome_number ; $index1++ ) {
	for ( my $index2 = 0 ; $index2 < $genome_number ; $index2++ ) {
	    $pairwise_adj_breaks[$index1][$index2] = 0;
	}
    }

    print STDERR "Calculating adjacency matrix for level $level ...\n";
    print STDERR "Printing adjacency vector for level $level ...\n";
    open (my $corevecfile, ">$outprefix" . "$level" . "_core_adjacency_vector.txt");
    open (my $noncorevecfile, ">$outprefix" . "$level" . "_noncore_adjacency_vector.txt");
    my $breaksvectorfile;
    if ($level == 100) {
	open ($breaksvectorfile, ">$outprefix" . "$level" . "_core_adjacency_breaks_vector.txt");
    }
    foreach my $cluster_id (1 .. $#cluster_size)  {
	my %adjacency_vector = (); # Key1 = cluster_number_1+end(_5 or _3), Key2 = cluster_number_2+end(_5 or_3), Value = array of genomes, 1 if the edge is present in the genome
	my $cluster_id_5 = $cluster_id . '_5';
	my $cluster_id_3 = $cluster_id . '_3';
	foreach my $cluster_member ( @{ $cluster_for_number[$cluster_id] } ) {
	    my $feat_id = $cluster_member->{'id'};
	    my $feat_orient = $feat_hash{$feat_id}->{'orient'};
	    my $next_3 = undef;
	    my $next_5 = undef;
	    my $feat_tag = $FeatnameLookupTag_hash{$feat_id};
	    my $genome_index = $TagIndex{$feat_tag};
	    my $feat_asmbl_id = $AssemblyLookup_hash{$feat_id};
	    my $AssemblyIndex = $TagByPointer{$feat_id}; # points to the position in the genome_hash array
	    my $i;
	    my $AssemblyIndexLast = $#{ $genome_hash{$feat_tag}{$feat_asmbl_id} };
	    if (!defined $AssemblyIndex) {
		print STDERR "ERROR: bad cluster $feat_id $feat_orient $feat_tag $feat_asmbl_id $AssemblyIndex $AssemblyIndexLast\n";
		next;
	    }
	    if ($DEBUG) {
		print STDERR "$feat_id $feat_orient $feat_tag $feat_asmbl_id $AssemblyIndex $AssemblyIndexLast\n";
	    }
	    for ($i = $AssemblyIndex + 1; $i <= $AssemblyIndexLast; $i++) { # iterate upward over adjacent proteins
		if (defined $genome_hash_context{$feat_tag}{$feat_asmbl_id}->{$i}) {
		    last; # stop at CONTEXT nodes
		}
		my $tmp_feat_id = ${ $genome_hash{$feat_tag}{$feat_asmbl_id} }[$i];
		if ($cluster_size[$cluster_number{$tmp_feat_id}] >= $size_cutoff) {
		    $next_3 = $cluster_number{$tmp_feat_id};
		    if ($feat_hash{$tmp_feat_id}->{'orient'} == 1) {#protein is on forward strand
			$next_3 = $next_3 . '_5';
		    } else {
			$next_3 = $next_3 . '_3';
		    }
		    last;
		}
	    }   
	    for ($i = $AssemblyIndex - 1; $i >= 0; $i--) { # iterate downward over adjacent proteins
		if (defined $genome_hash_context{$feat_tag}{$feat_asmbl_id}->{$i}) {
		    last; # stop at CONTEXT nodes
		}
		my $tmp_feat_id = ${ $genome_hash{$feat_tag}{$feat_asmbl_id} }[$i];
		if ($cluster_size[$cluster_number{$tmp_feat_id}] >= $size_cutoff) {
		    $next_5 = $cluster_number{$tmp_feat_id};
		    if ($feat_hash{$tmp_feat_id}->{'orient'} == 1) {#protein is on forward strand
			$next_5 = $next_5 . '_3';
		    } else {
			$next_5 = $next_5 . '_5';
		    }
		    last;
		}
	    }
	    if ($feat_orient == -1) {#protein is on reverse strand so swap ends
		my $tmp = $next_3;
		$next_3 = $next_5;
		$next_5 = $tmp;
	    }
	    if (defined $next_5) {
		if (defined $adjacency_matrix{$cluster_id_5}{$next_5}) {
		    if ($adjacency_vector{$cluster_id_5}{$next_5}->[$genome_index] == 0) {
			$adjacency_matrix{$cluster_id_5}{$next_5}++;
			$adjacency_vector{$cluster_id_5}{$next_5}->[$genome_index] = 1;
		    } else {#if not 0 then we've already seen this edge due to CONTEXT nodes
			die ("ERROR: should have avoided CONTEXT nodes above $cluster_id_5 $next_5\n");
		    }
		} else {
		    $adjacency_matrix{$cluster_id_5}{$next_5} = 1;
		    $adjacency_vector{$cluster_id_5}{$next_5} = [];
		    for ( my $index = 0 ; $index < $genome_number ; $index++ ) {
			$adjacency_vector{$cluster_id_5}{$next_5}->[$index] = 0;
		    }
		    $adjacency_vector{$cluster_id_5}{$next_5}->[$genome_index] = 1;
		}
	    }
	    if (defined $next_3) {
		if (defined $adjacency_matrix{$cluster_id_3}{$next_3}) {
		    if ($adjacency_vector{$cluster_id_3}{$next_3}->[$genome_index] == 0) {
			$adjacency_matrix{$cluster_id_3}{$next_3}++;
			$adjacency_vector{$cluster_id_3}{$next_3}->[$genome_index] = 1;
		    } else {#if not 0 then we've already seen this edge due to CONTEXT nodes
			die ("ERROR: should have avoided CONTEXT nodes above $cluster_id_3 $next_3\n");
		    }
		} else {
		    $adjacency_matrix{$cluster_id_3}{$next_3} = 1;
		    $adjacency_vector{$cluster_id_3}{$next_3} = [];
		    for ( my $index = 0 ; $index < $genome_number ; $index++ ) {
			$adjacency_vector{$cluster_id_3}{$next_3}->[$index] = 0;
		    }
		    $adjacency_vector{$cluster_id_3}{$next_3}->[$genome_index] = 1;
		}
	    }
	}
	my $diff_break = 1;
	for ( my $index1 = 0 ; $index1 < $genome_number ; $index1++ ) {
	    $vector_breaks[$index1] = "";
	}
	foreach my $next_5 (keys %{ $adjacency_vector{$cluster_id_5} } )  {
	    if ($cluster_size[$cluster_id] >= $size_cutoff) {
		print $corevecfile "($cluster_id_5,$next_5)";
		for ( my $index1 = 0 ; $index1 < $genome_number ; $index1++ ) {
		    print $corevecfile "\t$adjacency_vector{$cluster_id_5}{$next_5}->[$index1]";
		    if ($adjacency_vector{$cluster_id_5}{$next_5}->[$index1] > 0) {
			$vector_breaks[$index1] = $next_5;
		    }
		}
		print $corevecfile "\n";
		$diff_break++;
	    } else {
		print $noncorevecfile "($cluster_id_5,$next_5)";
		for ( my $index = 0 ; $index < $genome_number ; $index++ ) {
		    print $noncorevecfile "\t$adjacency_vector{$cluster_id_5}{$next_5}->[$index]";
		}
		print $noncorevecfile "\n";
	    }
	}
	if ($level == 100) {
	    $diff_break -= 2;
	    if ($diff_break > 0) {
		$hash_breaks{$cluster_id_5} = $diff_break;
		print $breaksvectorfile "$cluster_id_5";
		for ( my $index1 = 0 ; $index1 < $genome_number ; $index1++ ) {
		    print $breaksvectorfile "\t$vector_breaks[$index1]";
		}
		print $breaksvectorfile "\n";
	    }
	}
	for ( my $index1 = 0 ; $index1 < $genome_number ; $index1++ ) {
	    for ( my $index2 = 0 ; $index2 < $genome_number ; $index2++ ) {
		#if (($vector_breaks[$index1] ne "") && ($vector_breaks[$index2] ne "") && ($vector_breaks[$index1] ne $vector_breaks[$index2])) {
		    # do not count missing adjacencies where one or both of the entries is null
		#unfortunately the line above was not a distance metric (failed the triangular inequality for starts) was trying not to penalize missing data
		if ($vector_breaks[$index1] ne $vector_breaks[$index2]) {
		    $pairwise_adj_breaks[$index1][$index2]++;
		    if ($DEBUG) {
			print STDERR "break $cluster_id_5 $index1,$index2 $vector_breaks[$index1]:$vector_breaks[$index2]\n";
		    }
		}
	    }
	}
	$diff_break = 1;
	for ( my $index1 = 0 ; $index1 < $genome_number ; $index1++ ) {
	    $vector_breaks[$index1] = "";
	}
	foreach my $next_3 (keys %{ $adjacency_vector{$cluster_id_3} } )  {
	    if ($cluster_size[$cluster_id] >= $size_cutoff) {
		print $corevecfile "($cluster_id_3,$next_3)";
		for ( my $index1 = 0 ; $index1 < $genome_number ; $index1++ ) {
		    print $corevecfile "\t$adjacency_vector{$cluster_id_3}{$next_3}->[$index1]";
		    if ($adjacency_vector{$cluster_id_3}{$next_3}->[$index1] > 0) {
			$vector_breaks[$index1] = $next_3;
		    }
		}
		print $corevecfile "\n";
		$diff_break++;
	    } else {
		print $noncorevecfile "($cluster_id_3,$next_3)";
		for ( my $index = 0 ; $index < $genome_number ; $index++ ) {
		    print $noncorevecfile "\t$adjacency_vector{$cluster_id_3}{$next_3}->[$index]";
		}
		print $noncorevecfile "\n";
	    }
	}
	if ($level == 100) {
	    $diff_break -= 2;
	    if ($diff_break > 0) {
		$hash_breaks{$cluster_id_3} = $diff_break;
		print $breaksvectorfile "$cluster_id_3";
		for ( my $index1 = 0 ; $index1 < $genome_number ; $index1++ ) {
		    print $breaksvectorfile "\t$vector_breaks[$index1]";
		}
		print $breaksvectorfile "\n";
	    }
	}
	for ( my $index1 = 0 ; $index1 < $genome_number ; $index1++ ) {
	    for ( my $index2 = 0 ; $index2 < $genome_number ; $index2++ ) {
		#if (($vector_breaks[$index1] ne "") && ($vector_breaks[$index2] ne "") && ($vector_breaks[$index1] ne $vector_breaks[$index2])) {
		    # do not count missing adjacencies where one or both of the entries is null
		#unfortunately the line above was not a distance metric (failed the triangular inequality for starts) was trying not to penalize missing data
		if ($vector_breaks[$index1] ne $vector_breaks[$index2]) {
		    $pairwise_adj_breaks[$index1][$index2]++;
		    if ($DEBUG) {
			print STDERR "break $cluster_id_3 $index1,$index2 $vector_breaks[$index1]:$vector_breaks[$index2]\n";
		    }
		}
	    }
	}
    }
    close ($corevecfile);
    close ($noncorevecfile);

    if ($level == 100) {
	open (my $breaksfile, ">$outprefix" . "$level" . "_core_adjacency_breaks.txt");
	my $sort_split_cluster_id = sub {#sort routine used directly below to split the cluster id and sort numerically ascending
	    my ($cluster_ida, $enda) = split ('_', $a, 2);
	    my ($cluster_idb, $endb) = split ('_', $b, 2);
	    if ($cluster_ida < $cluster_idb){ #sort ascending numerical cluster id
		return -1;
	    } elsif ($cluster_ida > $cluster_idb){
		return 1;
	    } else { #sort descending by end 5 or 3 numerical
		return ($endb <=> $enda);
	    }
	};

	foreach my $cluster_end (sort $sort_split_cluster_id keys %hash_breaks) {
	    print $breaksfile "$cluster_end\t$hash_breaks{$cluster_end}\n";
	}
	close ($breaksfile);
	close ($breaksvectorfile);
    }

    open (my $adjdistpfile, ">$outprefix" . "$level" . "_core_adjacency_distance_matrix_phylip.txt");
    &print_pairwise_matrix (\@pairwise_adj_breaks, $adjdistpfile, 0.5, 0, 0, 1); #dividing by two because each break should get counted twice above unless there is a contig break
    close ($adjdistpfile);

    open (my $adjdistfile, ">$outprefix" . "$level" . "_core_adjacency_distance_matrix.txt");
    &print_pairwise_matrix (\@pairwise_adj_breaks, $adjdistfile, 0.5, 0, 0, 0); #dividing by two because each break should get counted twice above unless there is a contig break
    close ($adjdistfile);

    open (my $coreadjfile, ">$outprefix" . "$level" . "_core_adjacency_matrix.txt");
    foreach my $cluster_id (1 .. $#cluster_size)  {
	if ($cluster_size[$cluster_id] < $size_cutoff) {
	    next; #this loop outputs the core adjacency matrix from core to core clusters based on $level parameter
	}
	if (defined $cluster_visited{$cluster_id}) {
	    next; #skip this cluster - we already visited it following a previous core path
	}
	my $cluster_member = ${ $cluster_for_number[$cluster_id] }[0];
	my $feat_id = $cluster_member->{'id'};
	my $feat_orient = $feat_hash{$feat_id}->{'orient'};
	my $cur_cluster = $cluster_id;
	my $cur_strand = "+";
	if ($feat_orient == -1) {#protein is on reverse strand
	    $cur_strand = "-";
	}
	while (defined $cur_cluster) {
	    if (defined $cluster_visited{$cur_cluster}) {
		last; #stop following this path - we already visited it following a previous core path
	    }
	    $cluster_visited{$cur_cluster} = 1;
	    $cluster_member = ${ $cluster_for_number[$cur_cluster] }[0];
	    $feat_id = $cluster_member->{'id'};
	    print $coreadjfile "$cur_cluster$cur_strand:$cluster_size[$cur_cluster]\t$feat_id\t$feat_hash{$feat_id}->{'header'}\n";
	    if ($cur_strand eq "+") {
		&print_adjacency_matrix($coreadjfile, $cur_cluster . '_5', \%adjacency_matrix);
		print $coreadjfile " ||| ";
		($cur_cluster, $cur_strand) = &print_adjacency_matrix($coreadjfile, $cur_cluster . '_3', \%adjacency_matrix);
	    } else {
		&print_adjacency_matrix($coreadjfile, $cur_cluster . '_3', \%adjacency_matrix);
		print $coreadjfile " ||| ";
		($cur_cluster, $cur_strand) = &print_adjacency_matrix($coreadjfile, $cur_cluster . '_5', \%adjacency_matrix);
	    }
	    print $coreadjfile "\n";
	}
    }
    close ($coreadjfile);

    open (my $noncoreadjfile, ">$outprefix" . "$level" . "_noncore_adjacency_matrix.txt");
    foreach my $cluster_id (1 .. $#cluster_size)  {
	if ($cluster_size[$cluster_id] >= $size_cutoff) {
	    next; #this loop outputs the noncore adjacency matrix from noncore to core clusters based on $level parameter
	}
	if (defined $cluster_visited{$cluster_id}) {
	    die ("ERROR $cluster_id : $cluster_size[$cluster_id] : noncore clusters should not be visited!\n");
	}
	my $cluster_member = ${ $cluster_for_number[$cluster_id] }[0];
	my $feat_id = $cluster_member->{'id'};
	my $feat_orient = $feat_hash{$feat_id}->{'orient'};
	my $cur_cluster = $cluster_id;
	my $cur_strand = "+";
	if ($feat_orient == -1) {#protein is on reverse strand
	    $cur_strand = "-";
	}
	print $noncoreadjfile "$cur_cluster$cur_strand:$cluster_size[$cur_cluster]\t$feat_id\t$feat_hash{$feat_id}->{'header'}\n";
	if ($cur_strand eq "+") {
	    &print_adjacency_matrix($noncoreadjfile, $cur_cluster . '_5', \%adjacency_matrix);
	    print $noncoreadjfile " ||| ";
	    &print_adjacency_matrix($noncoreadjfile, $cur_cluster . '_3', \%adjacency_matrix);
	} else {
	    &print_adjacency_matrix($noncoreadjfile, $cur_cluster . '_3', \%adjacency_matrix);
	    print $noncoreadjfile " ||| ";
	    &print_adjacency_matrix($noncoreadjfile, $cur_cluster . '_5', \%adjacency_matrix);
	}
	print $noncoreadjfile "\n";
    }
    close ($noncoreadjfile);

}

sub calc_adjacency { # construct adjacency matrix for each level given

    foreach my $level (@cluster_levels)  {
	if ($DEBUG) {
	    print STDERR "calc_adjacency for level $level\n";
	}
	&calc_adjacency_weights($level);
    }

}

sub print_pairwise_cluster_matches { # print the match characteristics for matches within clusters and between clusters

    
    open (OUTPAIRWISEIN, ">$outprefix" . "pairwise_in_cluster.txt");

    print OUTPAIRWISEIN "#clust\tqtag\tqid\tstag\tsid\t%id\tevalue\tbitscor\tnorm\tBSR\tbest\tbibest\tsynbest\tsynbibe\tCGNbibe\tfull\tanchor\textend\tcliquet\tcliquea\n";
    foreach my $cluster_id (1 .. $#cluster_size)  {
	my @shift_cluster_array = @{ $cluster_for_number[$cluster_id] };
	foreach my $cluster_member_1 ( @{ $cluster_for_number[$cluster_id] } ) {
	    my $qid = $cluster_member_1->{'id'};
	    my $qtag = $FeatnameLookupTag_hash{$qid};
	    shift(@shift_cluster_array);
	    foreach my $cluster_member_2 ( @shift_cluster_array ) {
		my $sid = $cluster_member_2->{'id'};
		my $stag = $FeatnameLookupTag_hash{$sid};
		my $percent_id;
		my $rh_ptr = $relationship_hash{$qid}{$sid};
		print OUTPAIRWISEIN "$cluster_id\t$qtag\t$qid\t$stag\t$sid\t";
		if (defined $rh_ptr) {
		    printf OUTPAIRWISEIN "%5.2f\t%e\t%d\t%5.2f\t%5.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $rh_ptr->{'id'}, $rh_ptr->{'eval'}, $rh_ptr->{'score'}, $rh_ptr->{'norm_score'}, $rh_ptr->{'BSR'}, $rh_ptr->{'best'}, $rh_ptr->{'bibest'}, $rh_ptr->{'synbest'}, $rh_ptr->{'synbibest'}, $rh_ptr->{'CGN_bibest'}, $rh_ptr->{'full'}, $rh_ptr->{'anchor'}, $rh_ptr->{'extend'}, $rh_ptr->{'clique_all'}, $rh_ptr->{'clique_top'};
		} else {
		    printf OUTPAIRWISEIN "%5.2f\t%e\t%d\t%5.2f\t%5.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		}
	    }
	}
    }

    close (OUTPAIRWISEIN);
    
    open (OUTPAIRWISEOUT, ">$outprefix" . "pairwise_out_cluster.txt");

    print OUTPAIRWISEOUT "#clust1\tclust2\tqtag\tqid\tstag\tsid\t%id\tevalue\tbitscor\tnorm\tBSR\tbest\tbibest\tsynbest\tsynbibe\tCGNbibe\tfull\tanchor\textend\tcliquet\tcliquea\n";
    foreach my $qid (keys %relationship_hash)  { # go through all featnames
	my $qtag = $FeatnameLookupTag_hash{$qid};
	my $cluster_number_qid = $cluster_number{$qid};
	foreach my $sid (keys %{ $relationship_hash{$qid} } )  { # go through all relationships
	    if ($qid eq $sid) {
		next; #skip matches to same protein
	    }
	    if (!defined $cluster_number{$qid}) {
		die ("cluster number not defined for $qid\n");
	    }
	    if (!defined $cluster_number{$sid}) {
		die ("cluster number not defined for $sid\n");
	    }
	    my $stag = $FeatnameLookupTag_hash{$sid};
	    my $cluster_number_sid = $cluster_number{$sid};
	    if ($qtag eq $stag) {
		next; #skip matches to same genome
	    }
	    if ($cluster_number_qid == $cluster_number_sid) {
		next; # we already output these above
	    }
	    if ($qid gt $sid) {
		next; # only enter qid,sid pairs once not symmetrically since they indicate the same match
	    }
	    my $rh_ptr = $relationship_hash{$qid}{$sid};
	    print OUTPAIRWISEOUT "$cluster_number_qid\t$cluster_number_sid\t$qtag\t$qid\t$stag\t$sid\t";
	    printf OUTPAIRWISEOUT "%5.2f\t%e\t%d\t%5.2f\t%5.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $rh_ptr->{'id'}, $rh_ptr->{'eval'}, $rh_ptr->{'score'}, $rh_ptr->{'norm_score'}, $rh_ptr->{'BSR'}, $rh_ptr->{'best'}, $rh_ptr->{'bibest'}, $rh_ptr->{'synbest'}, $rh_ptr->{'synbibest'}, $rh_ptr->{'CGN_bibest'}, $rh_ptr->{'full'}, $rh_ptr->{'anchor'}, $rh_ptr->{'extend'}, $rh_ptr->{'clique_all'}, $rh_ptr->{'clique_top'};
	}
    }

    close (OUTPAIRWISEOUT);
}

sub calc_BSR_distances { # determine the pairwise distances between genomes based on BSR

#need to run through twice to first calculate average percent identity and mean max BSR before calculating cluster membership measure because
#a match could be missing between cluster members for two nonoverlapping fragments or matches below cutoffs which we do not want to include as a zero

    my ($level) = @_;

    my $size_cutoff = ($level * $genome_number) / 100;
    my @pairwise_num_clusters = ();# Index1 = genome tag index, Index2 = genome tag index, Value = number of clusters shared between the two genomes (later normalized by num_clusters)
    my @num_clusters = ();      # Index = genome tag index, Value = number of clusters in the genome

    foreach my $index1 (0 .. $#tag_array) {
	$num_clusters[$index1] = 0;
	foreach my $index2 (0 .. $#tag_array) {
	    $pairwise_num_clusters[$index1][$index2] = 0;
	    $ave_per_id[$index1][$index2] = 0;
	    $mean_max_BSR[$index1][$index2] = 0;
	}
    }

    foreach my $cluster_id (1 .. $#cluster_size)  {
	print STDERR "cluster $cluster_id\n" if ($DEBUG);
	if ($cluster_size[$cluster_id] < $size_cutoff){
	    next;
	}
	foreach my $cluster_member_1 ( @{ $cluster_for_number[$cluster_id] } ) {
	    my $qid = $cluster_member_1->{'id'};
	    my $qTagIndex = $cluster_member_1->{'tag'};
	    print STDERR "$qid\n" if ($DEBUG);
	    $num_clusters[$qTagIndex]++;
	    foreach my $cluster_member_2 ( @{ $cluster_for_number[$cluster_id] } ) {
		my $sid = $cluster_member_2->{'id'};
		my $sTagIndex = $cluster_member_2->{'tag'};
		print STDERR "\t$sid\t" if ($DEBUG);
		if (defined $relationship_hash{$qid}{$sid}) {
		    $pairwise_num_clusters[$qTagIndex][$sTagIndex]++;
		    $ave_per_id[$qTagIndex][$sTagIndex] += $relationship_hash{$qid}{$sid}->{'id'};
		    my $max_BSR = $relationship_hash{$qid}{$sid}->{'norm_score'};
		    if ($relationship_hash{$sid}{$qid}->{'norm_score'} > $max_BSR) {
			$max_BSR = $relationship_hash{$sid}{$qid}->{'norm_score'};
		    }
		    $mean_max_BSR[$qTagIndex][$sTagIndex] += $max_BSR;
		}
	    }
	}
    }

    foreach my $index1 (0 .. $#tag_array) {
	foreach my $index2 (0 .. $#tag_array) {
	    if ($pairwise_num_clusters[$index1][$index2] > 0) {
		$ave_per_id[$index1][$index2] /= $pairwise_num_clusters[$index1][$index2];
		$mean_max_BSR[$index1][$index2] /= $pairwise_num_clusters[$index1][$index2];
	    }
	}
    }
    open (my $pidentfile, ">$outprefix" . "$level" . "_pairwise_identity_matrix.txt");
    &print_pairwise_matrix (\@ave_per_id, $pidentfile, 1, 0, 0, 0);
    close ($pidentfile);

    open (my $pBSRfile, ">$outprefix" . "$level" . "_pairwise_BSR_matrix.txt");
    &print_pairwise_matrix (\@mean_max_BSR, $pBSRfile, 100, 0, 0, 0);
    close ($pBSRfile);

    open (my $pBSRdistpfile, ">$outprefix" . "$level" . "_pairwise_BSR_distance_matrix_phylip.txt");
    &print_pairwise_matrix (\@mean_max_BSR, $pBSRdistpfile, -100, 100, 0, 1);
    close ($pBSRdistpfile);

    open (my $pBSRdistfile, ">$outprefix" . "$level" . "_pairwise_BSR_distance_matrix.txt");
    &print_pairwise_matrix (\@mean_max_BSR, $pBSRdistfile, -100, 100, 0, 0);
    close ($pBSRdistfile);

}

sub calc_all_BSR_distances { # determine the pairwise distances between genomes based on cluster membership

    foreach my $level (@BSR_cluster_levels)  {
	if ($DEBUG) {
	    print STDERR "calc_BSR_distance for level $level\n";
	}
	&calc_BSR_distances($level);
    }
}

sub calc_pairwise_cluster_distances { # determine the pairwise distances between genomes based on cluster membership

    my @pairwise_num_clusters = ();# Index1 = genome tag index, Index2 = genome tag index, Value = number of clusters shared between the two genomes (intersection)
    my @num_clusters = ();         # Index = genome tag index, Value = number of clusters in the genome
    my @Sorenson = ();             # Same as pairwise_num_clusters but using Sorenson normalization
    my @Jaccard = ();              # Same as pairwise_num_clusters but using Jaccard normalization
    
    foreach my $index1 (0 .. $#tag_array) {
	$num_clusters[$index1] = 0;
	foreach my $index2 (0 .. $#tag_array) {
	    $pairwise_num_clusters[$index1][$index2] = 0;
	}
    }

    foreach my $cluster_id (1 .. $#cluster_size)  {
	print STDERR "cluster $cluster_id\n" if ($DEBUG);
	foreach my $cluster_member_1 ( @{ $cluster_for_number[$cluster_id] } ) {
	    my $qid = $cluster_member_1->{'id'};
	    my $qTagIndex = $cluster_member_1->{'tag'};
	    print STDERR "$qid\n" if ($DEBUG);
	    $num_clusters[$qTagIndex]++;
	    foreach my $cluster_member_2 ( @{ $cluster_for_number[$cluster_id] } ) {
		my $sid = $cluster_member_2->{'id'};
		my $sTagIndex = $cluster_member_2->{'tag'};
		print STDERR "\t$sid\t" if ($DEBUG);
		$pairwise_num_clusters[$qTagIndex][$sTagIndex]++;
	    }
	}
    }

    foreach my $index1 (0 .. $#tag_array) {
	if ($num_clusters[$index1] == 0) {
	    $num_clusters[$index1] = 1; # this shouldn't happen but just in case avoid division by zero later
	}
    }
    foreach my $index1 (0 .. $#tag_array) {
	foreach my $index2 (0 .. $#tag_array) {
	    #$pairwise_num_clusters[$index1][$index2] is the intersection of clusters for the two genomes
	    $Sorenson[$index1][$index2] = (200 * $pairwise_num_clusters[$index1][$index2]) / ($num_clusters[$index1] + $num_clusters[$index2]); # scale to be like percent identity and normalize by average number of clusters in both genomes - Sorenson similarity
	    $Jaccard[$index1][$index2] = (100 * $pairwise_num_clusters[$index1][$index2]) / (($num_clusters[$index1] + $num_clusters[$index2]) - $pairwise_num_clusters[$index1][$index2]); # scale to be like percent identity and normalize by union of clusters in both genomes - Jaccard similarity
	}
    }

    open (my $Ssimfile, ">$outprefix" . "Sorenson_pairwise_cluster_similarity_matrix.txt");
    &print_pairwise_matrix (\@Sorenson, $Ssimfile, 1, 0, 0, 0);
    close ($Ssimfile);

    open (my $Sdistpfile, ">$outprefix" . "Sorenson_pairwise_cluster_distance_matrix_phylip.txt");
    &print_pairwise_matrix (\@Sorenson, $Sdistpfile, -1, 100, 0, 1);
    close ($Sdistpfile);

    open (my $Sdistfile, ">$outprefix" . "Sorenson_pairwise_cluster_distance_matrix.txt");
    &print_pairwise_matrix (\@Sorenson, $Sdistfile, -1, 100, 0, 0);
    close ($Sdistfile);

    open (my $Jsimfile, ">$outprefix" . "Jaccard_pairwise_cluster_similarity_matrix.txt");
    &print_pairwise_matrix (\@Jaccard, $Jsimfile, 1, 0, 0, 0);
    close ($Jsimfile);

    open (my $Jdistpfile, ">$outprefix" . "Jaccard_pairwise_cluster_distance_matrix_phylip.txt");
    &print_pairwise_matrix (\@Jaccard, $Jdistpfile, -1, 100, 0, 1);
    close ($Jdistpfile);

    open (my $Jdistfile, ">$outprefix" . "Jaccard_pairwise_cluster_distance_matrix.txt");
    &print_pairwise_matrix (\@Jaccard, $Jdistfile, -1, 100, 0, 0);
    close ($Jdistfile);

}

sub calc_cluster_centroids { # determine the centroid of each cluster based on BSR scores

    
    open (my $centroidfile, ">$outprefix" . "centroids.fasta");
    foreach my $cluster_id (1 .. $#cluster_size)  {
	my $max_BSR_sum = -1;
	my $centroid = undef;
	print STDERR "cluster $cluster_id\n" if ($DEBUG);
	foreach my $cluster_member ( @{ $cluster_for_number[$cluster_id] } ) {
	    my $qid = $cluster_member->{'id'};
	    my $BSR_sum = 0;
	    print STDERR "$qid\n" if ($DEBUG);
	    foreach my $sid (keys %{ $relationship_hash{$qid} } )  {
		print STDERR "\t$sid\t" if ($DEBUG);
		if ($cluster_id != $cluster_number{$sid}) {
		    print STDERR "in $cluster_number{$sid} $relationship_hash{$qid}{$sid}->{'BSR'}\n" if ($DEBUG);
		    next; #skip matches outside of the cluster
		}
		print STDERR "$relationship_hash{$qid}{$sid}->{'BSR'}\n" if ($DEBUG);
		$BSR_sum += $relationship_hash{$qid}{$sid}->{'BSR'};
	    }
	    print STDERR "$BSR_sum\n" if ($DEBUG);
	    if ($BSR_sum > $max_BSR_sum) {
		$max_BSR_sum = $BSR_sum;
		$centroid = $qid;
	    }
	}
	print STDERR "$max_BSR_sum\n" if ($DEBUG);
	if (!defined $centroid) {
	    die "centroid could not be determined for cluster $cluster_id !\n";
	}
	$centroids[$cluster_id] = $centroid;
	print $centroidfile ">centroid_$cluster_id $centroid $feat_hash{$centroid}->{'header'}\n";
	my $seq_len = $feat_hash{$centroid}->{'length'};
	my $sequence = $feat_hash{$centroid}->{'sequence'};
	for ( my $pos = 0 ; $pos < $seq_len ; $pos += 60 ) {
	    print $centroidfile substr($sequence, $pos, 60), "\n";
	}

    }
    close ($centroidfile);

}

sub write_files  {

    my ($choice, $query_column, $subject_column, $query_featname, $subject_tag, $subject_featname, $cluster_number) = @_;
    my $ratio = "";
    my $x = "";

    if ($print_cluster_numbers && ($subject_column == 1)) {
	    print OUTVENN "$cluster_number\t" if ($vennfile);
	    print OUTZEROONE "$cluster_number\t" if ($vennfile);
	    print OUTVENNID "$cluster_number\t" if ($vennfile);
	    if ($query_column == 1)  {
		printf MICRO "$cluster_number\t" if ($microarray);
		print PAULSEN "$cluster_number\t" if ($hitsfile);
		print READ "$cluster_number\t" if ($normalizefile);
		printf OUTID "$cluster_number\t" if ($vennfile);
	    }
    }
    if ($choice == "1")  {
	if ($query_column == 1)  {
	    print OUTVENN "$query_featname" if ($vennfile);
	    print MICRO "$query_featname\t$feat_hash{$query_featname}->{'header'}" if ($microarray);
	    print PAULSEN "$query_featname \[$feat_hash{$query_featname}->{'header'}\]" if ($hitsfile);
	    print READ "$query_featname\t$feat_hash{$query_featname}->{'header'}" if ($normalizefile);
	    print OUTZEROONE "1" if ($vennfile);
	    print OUTVENNID "$query_featname" if ($vennfile);
	    print OUTID "$query_featname\t$feat_hash{$query_featname}->{'header'}" if ($vennfile);
	    
	}
	else  {
	    print OUTVENN "\t$query_featname" if ($vennfile);
	    print OUTZEROONE "\t1" if ($vennfile);
	    print OUTVENNID "\t$query_featname" if ($vennfile);
	}
    } elsif ($choice == "2")  {
	my $score;
	my $percent_id;
	if (!defined $relationship_hash{$query_featname}{$subject_featname}) {
	    print STDERR "undefined relationship_hash for $query_featname and $subject_featname in the same cluster\n" if ($DEBUG);
	    $score = 0;
	    $percent_id = 0;
	} else {
	    $score = $relationship_hash{$query_featname}{$subject_featname}->{'score'};
	    $percent_id = $relationship_hash{$query_featname}{$subject_featname}->{'id'};
	}
        if ($subject_column != 1)  {
	    print OUTVENN "\t" if ($vennfile);
	    print OUTZEROONE "\t" if ($vennfile);
	    print OUTVENNID "\t" if ($vennfile);
	    if ($query_column == 1) { # only print the reference as query, not the leftovers from the other genomes
		printf MICRO "\t" if ($microarray);
		print PAULSEN "\t" if ($hitsfile);
		print READ "\t" if ($normalizefile);
		printf OUTID "\t" if ($vennfile);
	    }
	}
	print OUTVENN "$subject_featname" if ($vennfile);
	print OUTZEROONE "1" if ($vennfile);
	printf OUTVENNID "$subject_featname (%5.2f", $percent_id if ($vennfile);
	print OUTVENNID "%)" if ($vennfile);
	if ($query_column == 1) { # only print the reference as query, not the leftovers from the other genomes
            #$Qbytaghash{$qid}{$stag}->{$sid}->{'BSR'}
	    $x = $score/$relationship_hash{$query_featname}{$query_featname}->{'score'};  # determine normalized BLAST score
	    ## perhaps this can be simplified by getting data from Qbytaghash?
	    $ratio = -99*$x+100; # convert to scale 1 = perfect match, 100 = no match ( arthematic supplied by Emmanuel Mongodin )
	    printf MICRO "%2.2f", $ratio if ($microarray);
	    print PAULSEN "$subject_featname \[$feat_hash{$subject_featname}->{'header'}\]" if ($hitsfile);
	    print READ "$x" if ($normalizefile);
	    printf OUTID "%5.2f", $percent_id if ($vennfile);
	}

    } elsif ($choice == "3")  {
	if ($subject_column != 1)  {
	    print OUTVENN "\t" if ($vennfile); 
	    print OUTZEROONE "\t" if ($vennfile);
	    print OUTVENNID "\t" if ($vennfile);
	    if ($query_column == 1) { # only print the reference as query, not the leftovers from the other genomes
                printf MICRO "\t" if ($microarray);
                print PAULSEN "\t" if ($hitsfile);
                print READ "\t" if ($normalizefile);
                print OUTID "\t" if ($vennfile);
	    }
	}
	print OUTZEROONE "0" if ($vennfile);
	print OUTVENN "----------" if ($vennfile);  # if no match print the default delimiter ----------
	if ($query_column == 1) { # only print the reference as query, not the leftovers from the other genomes
	    $x = 0;
	    $ratio = -99*$x+100;
	    printf MICRO "%2.2f", $ratio if ($microarray);
	    print READ "$x" if ($normalizefile);
	}
    } elsif ($choice == "4")  {
	print OUTVENN "\n" if ($vennfile);
	print OUTZEROONE "\n" if ($vennfile);
	print OUTVENNID "\n" if ($vennfile);
	if ($query_column == 1)  {
	    print MICRO "\n" if ($microarray);
	    print PAULSEN "\n" if ($hitsfile);
	    print READ "\n" if ($normalizefile);
	    print OUTID "\n" if ($vennfile);
	}
    }
}

sub synteny {#&synteny($genome_tag, $query_featname, $subject_tag, $subject_featname)
# need to make genome_hash two dimensional to account for multiple genomes and multiple molecules
    my ($query_tag, $query_featname, $subject_tag, $subject_featname, $window) = @_;
    my $QueryArrayIndex = $TagByPointer{$query_featname}; # points to the position in the genome_hash array
    my $SubjectArrayIndex = $TagByPointer{$subject_featname}; # points to the position in the genome_hash array
    my $query_orient = $feat_hash{$query_featname}->{'orient'};
    my $subject_orient = $feat_hash{$subject_featname}->{'orient'};
    my $orient = $query_orient * $subject_orient;
    my $i = "";
    my $j = "";
    my $OuterLoopEnd = "";
    my $OuterLoopBeg = "";
    my $InnerLoopEnd = "";
    my $InnerLoopBeg = "";
    my $query = "";
    my $subject = "";
    my $Qasmbl_id = $AssemblyLookup_hash{$query_featname};
    my $Sasmbl_id = $AssemblyLookup_hash{$subject_featname};    
    my $query_array_last = $#{ $genome_hash{$query_tag}{$Qasmbl_id} };
    my $subject_array_last = $#{ $genome_hash{$subject_tag}{$Sasmbl_id} };
    my $synteny_range = $window; #this is the number of proteins on either side of the query and subject to compute synteny
    my $num_synteny = 0;
    my $num_bibest_window = 0;
    my $individual_score;
    my $max_individual_score;
    my $subject_bibest;
    my $max_bibest;
    my $total_score = 0;
    my %max_subject_hash = ();
    my $max_subject;
    my $max_j;
    my $abs_rel_pos_query;
    my $abs_rel_pos_subject;
    my $dist_from_center;
    my $query_5pend; #+1 if on the 5' end side of the subject ORF -1 3'end
    my $subject_5pend; #+1 if on the 5' end side of the subject ORF -1 3'end


    print STDERR "query $query_featname ($QueryArrayIndex)     subject $subject_featname ($SubjectArrayIndex)\n" if ($DEBUG);

    $OuterLoopEnd = $QueryArrayIndex + $synteny_range;
    if ($OuterLoopEnd > $query_array_last) {# check to make sure we don't extend past the last array element 
	$OuterLoopEnd = $query_array_last; # if so, make OuterLoopEnd equal to the size of the array
    }
    $OuterLoopBeg = $QueryArrayIndex - $synteny_range;
    if ($OuterLoopBeg < 0)  { # don't let this run off the beginning of the array
	$OuterLoopBeg = 0; # adjust to the first position of the array
    }
    for ($i = $OuterLoopBeg; $i <= $OuterLoopEnd; $i++) { # iterate over synteny range of query genome
	print STDERR "       qpos = $i (${ $genome_hash{$query_tag}{$Qasmbl_id} }[$i])\n" if ($DEBUG);
	$abs_rel_pos_query = abs ($i - $QueryArrayIndex);
	$query = ${ $genome_hash{$query_tag}{$Qasmbl_id} }[$i];
	if (($QueryArrayIndex - $i) * $query_orient > 0) {
	    $query_5pend = 1;
	} else {
	    $query_5pend = -1;
	}
	$InnerLoopEnd = $SubjectArrayIndex + $synteny_range;
	if ($InnerLoopEnd > $subject_array_last) {# check to make sure we don't extend past the last array element 
	    $InnerLoopEnd = $subject_array_last; # if so, make InnerLoopEnd equal to the size of the array
	}
	$InnerLoopBeg = $SubjectArrayIndex - $synteny_range;
	if ($InnerLoopBeg < 0)  { # don't let this run off the beginning of the array
	    $InnerLoopBeg = 0; # adjust to the first position of the array
	}
	$max_individual_score = 0;
	$max_bibest = 0;
	for ($j = $InnerLoopBeg; $j <= $InnerLoopEnd; $j++) { # iterate over synteny range of subject genome
	    $abs_rel_pos_subject = abs ($j - $SubjectArrayIndex);
	    if ($abs_rel_pos_query > $abs_rel_pos_subject) {
		$dist_from_center = $abs_rel_pos_query;
	    } else {
		$dist_from_center = $abs_rel_pos_subject;
	    }
	    $subject = ${ $genome_hash{$subject_tag}{$Sasmbl_id} }[$j];
	    if (($SubjectArrayIndex - $j) * $subject_orient > 0) {
		$subject_5pend = 1;
	    } else {
		$subject_5pend = -1;
	    }
	    $subject_bibest = 0;
	    if ((defined $relationship_hash{$query}{$subject}) && (defined $relationship_hash{$subject}{$query}) && (defined $relationship_hash{$query}{$subject}->{'BSR'}) && (defined $relationship_hash{$subject}{$query}->{'BSR'})) {#need to check for relationship_hash existence first to avoid creating it when we check for BSR existence
		$individual_score = 1;
		if ($relationship_hash{$query}{$subject}->{'best'}) {
		    $individual_score += 2;
		    if ($relationship_hash{$query}{$subject}->{'bibest'}) {
			$subject_bibest = 1;
			$individual_score += 5;
			if ($relationship_hash{$query}{$subject}->{'clique_top'}) {
			    $individual_score += (5 * $relationship_hash{$query}{$subject}->{'clique_top'}) / $genome_number;
			    if ($relationship_hash{$query}{$subject}->{'clique_all'}) {
				$individual_score += (5 * $relationship_hash{$query}{$subject}->{'clique_all'}) / $genome_number;
			    }
			}
		    }
		}
		if ($relationship_hash{$query}{$subject}->{'extend'}) {
		    $individual_score += 10;
		    $subject_bibest = 1; # count these as bibest matches also
		}
		if ($relationship_hash{$query}{$subject}->{'anchor'}) {
		    $individual_score += 20;
		    $subject_bibest = 1; # count these as bibest matches also
		}
		if (($orient * ($feat_hash{$query}->{'orient'} * $feat_hash{$subject}->{'orient'})) > 0) {
		    $individual_score += 1;
		}
		if ($query_5pend * $subject_5pend > 0) {
		    $individual_score += 1;
		}
		$individual_score += (($synteny_range + 1) - abs ($abs_rel_pos_query - $abs_rel_pos_subject));
		$individual_score *= exp ((log 3) *(($synteny_range + 1) - $dist_from_center));
		if ($individual_score > $max_individual_score) {
		    $max_individual_score = $individual_score;
		    $max_bibest = $subject_bibest;
		    $max_subject = $subject;
		    $max_j = $j;
		}
		print STDERR "         spos = $j (${ $genome_hash{$subject_tag}{$Sasmbl_id} }[$j]) [$relationship_hash{$query}{$subject}->{'id'}] <=> [$relationship_hash{$subject}{$query}->{'id'}] $individual_score $relationship_hash{$query}{$subject}->{'best'} $relationship_hash{$query}{$subject}->{'bibest'} $relationship_hash{$query}{$subject}->{'clique_top'} $relationship_hash{$query}{$subject}->{'clique_all'} (orient: $orient $feat_hash{$query}->{'orient'} $feat_hash{$subject}->{'orient'}) (5p: $query_5pend $subject_5pend)\n" if ($DEBUG);
	    }
	}
	$total_score += $max_individual_score;
	if ($max_individual_score > 0) {
	    if ($max_bibest) { # keep track of how many bibest matches are in CGN window including match being evaluated
		$num_bibest_window++;
	    }
	    if (!defined $max_subject_hash{$max_subject}) {
		$max_subject_hash{$max_subject} = $max_individual_score;
		if (($i != $QueryArrayIndex) && ($max_j != $SubjectArrayIndex)) { # don't count synteny for the match being evaluated
		    $num_synteny++;
		}
	    } else {
		if ($max_subject_hash{$max_subject} >= $max_individual_score) {
		    $total_score -= $max_individual_score;
		} else {
		    $total_score -= $max_subject_hash{$max_subject};
		    $max_subject_hash{$max_subject} = $max_individual_score;
		}
	    }
	}
    }
    print STDERR "     num_synteny $num_synteny num_bibest_window $num_bibest_window\n" if ($DEBUG);
    print STDERR "     total_score $total_score\n" if ($DEBUG);
    return($total_score, $num_bibest_window);

}

#GGS new subroutine to detect frameshifts/protein fragments and record them
#need to consider only doing this based on bidirectionally best blast matches
#would be better to take orientation and order of matches against the query into account
sub frameshifts {

    my %max_frames = ();        #value is counter for best/longest fragment of a frameshifted set, Key = feat_name
    my %con_comp_frames = ();   #value is counter for sets of frameshifted fragments, Keys are feat_name, feat_name

    my $qid = ""; # query feat_name
    my $sid = ""; # subject feat_name
    my $qtag = ""; # query genome tag
    my $stag = ""; # subject genome tag

    foreach $qtag (@tag_array)  {  # start looping through genomes by order in tag file
	foreach $qid (keys %{ $Tagbyfeatnamehash{$qtag} } ) { # go through featnames of each genome to look for matches in other genomes
	    foreach $stag (keys %{ $Qbytaghash{$qid} } ) {# if query protein matches anything in subject genome, lets drill through each match

		my $sort_asmbl_id_genome_order = sub {#sort routine used in frameshifts to sort proteins from a given genome by asmbl_id and then by position in the assembly/contig
		    my $asmbl_id_a = $AssemblyLookup_hash{$a};
		    my $asmbl_id_b = $AssemblyLookup_hash{$b};
		    if ($asmbl_id_a lt $asmbl_id_b){ #sort ascending asmbl_id
			return -1;
		    } elsif ($asmbl_id_a gt $asmbl_id_b){
			return 1;
		    } else { #sort ascending genome order
			$TagByPointer{$a} <=> $TagByPointer{$b};
		    }
		};

		my @qid_array_of_sids = sort $sort_asmbl_id_genome_order keys %{ $Qbytaghash{$qid}{$stag} };
		my $cur_index = 0;
		my $skip_index = 0;

		foreach my $sid (@qid_array_of_sids) {
		    if (($qid eq $sid) || ($relationship_hash{$qid}{$sid}->{"id"} < $fspercentid)) { # remove low %ID matches and self match
			$skip_index++;
		    } else {
			$qid_array_of_sids[$cur_index] = $qid_array_of_sids[$skip_index];
			$cur_index++;
			$skip_index++;
		    }
		}
		$cur_index--;
		$#qid_array_of_sids = $cur_index;
		$cur_index = 0;

		if (@qid_array_of_sids < 2) {
		    next;
		}# only need to look for frameshifts with two or more matches
		print STDERR "frameshifts: $qtag $qid $stag\n" if ($DEBUG);
		print STDERR "@qid_array_of_sids\n" if ($DEBUG);

		my $prev_asmbl_id = "";
		my $cur_asmbl_id;
		my $query_hits_span;
		my $min_query;
		my $max_query;
		my $cur_min_query;
		my $cur_max_query;
		my $cur_query_hit_len;
		my $sum_query_hit_len;
		my $start_index = 0;
		my $stop_index = 0;
		my $max_index;
		my $max_query_hit_len;
		my $failed_to_extend;
		my $prev_pos;
		my $cur_pos;

		if (@qid_array_of_sids == 2) {#look for fragments near the ends of two contigs (could be the same contig for circular contigs)
		    my $sid0 = $qid_array_of_sids[0];
		    my $sid1 = $qid_array_of_sids[1];
		    my $asmbl_id0 = $AssemblyLookup_hash{$sid0};
		    my $asmbl_id1 = $AssemblyLookup_hash{$sid1};
		    my $pos0 = $TagByPointer{$sid0};
		    my $pos1 = $TagByPointer{$sid1};
		    my $last_pos0 = $#{ $genome_hash{$stag}{$asmbl_id0} };
		    my $last_pos1 = $#{ $genome_hash{$stag}{$asmbl_id1} };
		    if ((($pos0 <= 1) || ($pos0 >= ($last_pos0 - 1))) && (($pos1 <= 1) || ($pos1 >= ($last_pos1 - 1)))) {
			#matches are near ends of contigs so treat as possible fragments
			#if we've added CONTEXT for circular contigs/chromosomes this will not work
			print STDERR "matches are near ends of contigs\n" if ($DEBUG);

			$min_query = $relationship_hash{$qid}{$sid0}->{'min_query'};
			$max_query = $relationship_hash{$qid}{$sid0}->{'max_query'};
			$query_hits_span = ($max_query - $min_query) + 1;
			$sum_query_hit_len = $query_hits_span;
			$max_query_hit_len = $query_hits_span;
			$max_index = 0;
			
			$cur_min_query = $relationship_hash{$qid}{$sid1}->{'min_query'};
			if ($cur_min_query < $min_query) {
			    $min_query = $cur_min_query;
			}
			$cur_max_query = $relationship_hash{$qid}{$sid1}->{'max_query'};
			if ($cur_max_query > $max_query) {
			    $max_query = $cur_max_query;
			}
			$query_hits_span = ($max_query - $min_query) + 1;
			$cur_query_hit_len = ($cur_max_query - $cur_min_query) + 1;
			$sum_query_hit_len += $cur_query_hit_len;
			if ($sum_query_hit_len < ($frameshift_length_ratio * $query_hits_span)){
			    if ($cur_query_hit_len > $max_query_hit_len){
				$max_query_hit_len = $cur_query_hit_len;
				$max_index = 1;
			    }
			    #record frameshift/fragments
			    print STDERR "recording as frameshift\n" if ($DEBUG);
			    $feat_hash{$sid0}->{'fragment'} += 2; # full length matches are doubled so double these
			    $feat_hash{$sid1}->{'fragment'} += 2; # full length matches are doubled so double these
			    $feat_hash{$qid}->{'fusion'} += 2; # full length matches are doubled so double these
			    if (!defined $con_comp_frames{$sid0}{$sid1}){
				$con_comp_frames{$sid0}{$sid1} = 1;
				$con_comp_frames{$sid1}{$sid0} = 1;
			    } else {
				$con_comp_frames{$sid0}{$sid1}++;
				$con_comp_frames{$sid1}{$sid0}++;
			    }
			    if (!defined $max_frames{$sid0}){
				$max_frames{$sid0} = 0;
			    }
			    if (!defined $max_frames{$sid1}){
				$max_frames{$sid1} = 0;
			    }
			    if (!defined $max_frames{$qid_array_of_sids[$max_index]}){
				$max_frames{$qid_array_of_sids[$max_index]} = 1;
			    } else {
				$max_frames{$qid_array_of_sids[$max_index]}++;
			    }
			}
			next; #already tested for fragments/frameshifts if near the ends of contigs
		    }
		}

		foreach $sid (@qid_array_of_sids) {
		    $cur_asmbl_id = $AssemblyLookup_hash{$sid};
		    $cur_pos = $TagByPointer{$sid};

		    print STDERR "c: $cur_asmbl_id $sid $cur_pos p: $prev_asmbl_id\n" if ($DEBUG);

		    if ($cur_asmbl_id eq $prev_asmbl_id){
			$cur_min_query = $relationship_hash{$qid}{$sid}->{'min_query'};
			if ($cur_min_query < $min_query) {
			    $min_query = $cur_min_query;
			}
			$cur_max_query = $relationship_hash{$qid}{$sid}->{'max_query'};
			if ($cur_max_query > $max_query) {
			    $max_query = $cur_max_query;
			}
			if (($cur_min_query <= $max_missing_aa + 1) && ($cur_max_query + $max_missing_aa >= $feat_hash{$qid}->{'length'})) {
			    $failed_to_extend = 1; #not a protein fragment if almost a full length match
			#} elsif ((($prev_pos + 1) != $cur_pos) && (($prev_pos - 1) != $cur_pos)) {
			} elsif (abs($prev_pos - $cur_pos) > 3) {
			    $failed_to_extend = 1; #not protein fragments/frameshift if not near by on the contig (allow two intervening orfs for transposases)
			} else {
			    $query_hits_span = ($max_query - $min_query) + 1;
			    $cur_query_hit_len = ($cur_max_query - $cur_min_query) + 1;
			    $sum_query_hit_len += $cur_query_hit_len;
			    if ($sum_query_hit_len < ($frameshift_length_ratio * $query_hits_span)){
				$failed_to_extend = 0;
				if ($cur_query_hit_len > $max_query_hit_len){
				    $max_query_hit_len = $cur_query_hit_len;
				    $max_index = $cur_index;
				}
				$stop_index = $cur_index;
			    } else {
				$failed_to_extend = 1;
			    }
			}
			print STDERR "failed_to_exend = $failed_to_extend\n" if ($DEBUG);
		    }
		    
		    if (($cur_asmbl_id ne $prev_asmbl_id) || $failed_to_extend || ($cur_index == $#qid_array_of_sids)){
			if ($stop_index > $start_index) {#record the frameshifted featnames
			    for (; $start_index < $stop_index; $start_index++){
				print STDERR "recording as frameshift $qid_array_of_sids[$start_index] $qid_array_of_sids[$start_index + 1]\n" if ($DEBUG);
				$feat_hash{$qid_array_of_sids[$start_index]}->{'fragment'} += 2; # full length matches are doubled so double these
				$feat_hash{$qid_array_of_sids[$start_index + 1]}->{'fragment'} += 2; # full length matches are doubled so double these
				$feat_hash{$qid}->{'fusion'} += 2; # full length matches are doubled so double these
				if (!defined $con_comp_frames{$qid_array_of_sids[$start_index]}{$qid_array_of_sids[$start_index + 1]}){
				    $con_comp_frames{$qid_array_of_sids[$start_index]}{$qid_array_of_sids[$start_index + 1]} = 1;
				    $con_comp_frames{$qid_array_of_sids[$start_index + 1]}{$qid_array_of_sids[$start_index]} = 1;
				} else {
				    $con_comp_frames{$qid_array_of_sids[$start_index]}{$qid_array_of_sids[$start_index + 1]}++;
				    $con_comp_frames{$qid_array_of_sids[$start_index + 1]}{$qid_array_of_sids[$start_index]}++;
				}
				if (!defined $max_frames{$qid_array_of_sids[$start_index]}){
				    $max_frames{$qid_array_of_sids[$start_index]} = 0;
				}
			    }
			    if (!defined $max_frames{$qid_array_of_sids[$stop_index]}){
				$max_frames{$qid_array_of_sids[$stop_index]} = 0;
			    }
			    if (!defined $max_frames{$qid_array_of_sids[$max_index]}){
				$max_frames{$qid_array_of_sids[$max_index]} = 1;
			    } else {
				$max_frames{$qid_array_of_sids[$max_index]}++;
			    }
			}
			if ($cur_index == $#qid_array_of_sids){
			    last;
			}
			$min_query = $relationship_hash{$qid}{$sid}->{'min_query'};
			$max_query = $relationship_hash{$qid}{$sid}->{'max_query'};
			if (($min_query <= $max_missing_aa + 1) && ($max_query + $max_missing_aa >= $feat_hash{$qid}->{'length'})) {
			    $cur_asmbl_id = ""; #not a protein fragment if almost a full length match
			    print STDERR "Not a protein fragment: $min_query - $max_query; full length match: $feat_hash{$qid}->{'length'}\n" if ($DEBUG);
			} else {
			    $query_hits_span = ($max_query - $min_query) + 1;
			    $sum_query_hit_len = $query_hits_span;
			    $max_query_hit_len = $query_hits_span;
			    $max_index = $cur_index;
			    $start_index = $cur_index;
			    $stop_index = $cur_index;
			}
		    }

		    $prev_asmbl_id = $cur_asmbl_id;
		    $prev_pos = $cur_pos;
		    $cur_index++;
		}

	    }
	}
    }

    open (OUTFRAMESHIFT, ">$outprefix" . "frameshifts.txt");
    my $prev_asmbl_id = "";
    my $cur_asmbl_id;
    my $prev_tag = "";
    my $cur_tag;
    my $max_feat_id;
    my $sort_fragments = sub {#sort routine used in frameshifts to sort proteins genome tag, then asmbl_id and then by best/longest fragment
	my $taga = $FeatnameLookupTag_hash{$a};
	my $tagb = $FeatnameLookupTag_hash{$b};
	my $asmbl_id_a = $AssemblyLookup_hash{$a};
	my $asmbl_id_b = $AssemblyLookup_hash{$b};
	my $max_a = $max_frames{$a};
	my $max_b = $max_frames{$b};
	if ($taga lt $tagb){ #sort ascending genome tag
	    return -1;
	} elsif ($taga gt $tagb){
	    return 1;
	} elsif ($asmbl_id_a lt $asmbl_id_b){ #sort ascending asmbl_id
	    return -1;
	} elsif ($asmbl_id_a gt $asmbl_id_b){
	    return 1;
	} elsif ($max_a > $max_b) { #sort descending max_frames
	    return -1;
	} elsif ($max_a < $max_b) {
	    return 1;
	} else { #sort descending by length of protein
	    return ($feat_hash{$b}->{'length'} <=> $feat_hash{$a}->{'length'});
	}
    };

    foreach $max_feat_id (sort $sort_fragments keys %max_frames) { #sort protein fragments by genome tag, then asmbl_id and then by best/longest fragment
	if (!defined $max_frames{$max_feat_id}) {
	    print STDERR "already included max_feat_id = $max_feat_id in another set of fragments/frameshifts\n" if ($DEBUG);
	    next; #already included this in another set of fragments/frameshifts
	}
	print STDERR "max_feat_id = $max_feat_id : max_frames = $max_frames{$max_feat_id} :full $feat_hash{$max_feat_id}->{'full'}\n" if ($DEBUG);
	$cur_tag = $FeatnameLookupTag_hash{$max_feat_id};
	$cur_asmbl_id = $AssemblyLookup_hash{$max_feat_id};
	my @linked_array = ();
	foreach my $linked_feat_id (sort { $con_comp_frames{$max_feat_id}{$b} <=> $con_comp_frames{$max_feat_id}{$a} } keys %{ $con_comp_frames{$max_feat_id} }) {
	    #Add linked fragments that qualify to the linked_array
	    if ((!defined $max_frames{$linked_feat_id}) || (!defined $feat_hash{$linked_feat_id})) {
		delete $con_comp_frames{$max_feat_id}{$linked_feat_id};
		delete $con_comp_frames{$linked_feat_id}{$max_feat_id};
		next; #linked_feat_id was already deleted or marked as a fragment/frameshift
	    }
	    if (($con_comp_frames{$max_feat_id}{$linked_feat_id} < $frames_link_thresh) || ((2 * $con_comp_frames{$max_feat_id}{$linked_feat_id}) < $feat_hash{$max_feat_id}->{'full'}) || ((2 * $con_comp_frames{$max_feat_id}{$linked_feat_id}) < $feat_hash{$linked_feat_id}->{'full'})) {
		print STDERR "$linked_feat_id failed minimal linkage test\n" if ($DEBUG);
		next; #link less than threshold
	    }
	    push (@linked_array, $linked_feat_id);
	    #Remove links we have already used to add fragments
	    print STDERR "frameshift weight: $con_comp_frames{$max_feat_id}{$linked_feat_id}\n" if ($DEBUG);
	    print STDERR "linked_feat_id = $linked_feat_id :full $feat_hash{$linked_feat_id}->{'full'}\n" if ($DEBUG);
	    delete $con_comp_frames{$max_feat_id}{$linked_feat_id};
	    delete $con_comp_frames{$linked_feat_id}{$max_feat_id};
	}
	if (@linked_array < 1) {
	    print STDERR "all fragments failed minimal linkage test\n" if ($DEBUG);
	    next; #no valid fragments linked to $max_feat_id
	}
	if ($cur_tag ne $prev_tag){
	    print OUTFRAMESHIFT ">genome $cur_tag\n>asmbl_id $cur_asmbl_id\n";
	} elsif ($cur_asmbl_id ne $prev_asmbl_id) {
	    print OUTFRAMESHIFT ">asmbl_id $cur_asmbl_id\n";
	}
	print OUTFRAMESHIFT $max_feat_id;
	$feat_hash{$max_feat_id}->{'retained'} = 1;
	my $linked_feat_id;
	while (defined ($linked_feat_id = shift(@linked_array))) {
	    if (!defined $max_frames{$linked_feat_id}) {
		print STDERR "ERROR!!! already included $linked_feat_id in another set of fragments/frameshifts\n" if ($DEBUG);
		next; #already included this in another set of fragments/frameshifts
	    }
	    print STDERR "linked_feat_id = $linked_feat_id : max_frames = $max_frames{$linked_feat_id} :full $feat_hash{$linked_feat_id}->{'full'}\n" if ($DEBUG);
	    print OUTFRAMESHIFT "\t$linked_feat_id";
	    delete $max_frames{$linked_feat_id};

	    #Add all of the linked fragments to the @linked_aray for the $linked_feat_id
	    foreach my $new_linked_feat_id (sort { $con_comp_frames{$linked_feat_id}{$b} <=> $con_comp_frames{$linked_feat_id}{$a} } keys %{ $con_comp_frames{$linked_feat_id} }) {
		#Add linked fragments that qualify to the linked_array
		if ((!defined $max_frames{$new_linked_feat_id}) || (!defined $feat_hash{$new_linked_feat_id})) {
		    delete $con_comp_frames{$linked_feat_id}{$new_linked_feat_id};
		    delete $con_comp_frames{$new_linked_feat_id}{$linked_feat_id};
		    next; #new_linked_feat_id was already deleted or marked as a fragment/frameshift
		}
		if (($con_comp_frames{$linked_feat_id}{$new_linked_feat_id} < $frames_link_thresh) || ((2 * $con_comp_frames{$linked_feat_id}{$new_linked_feat_id}) < $feat_hash{$linked_feat_id}->{'full'}) || ((2 * $con_comp_frames{$linked_feat_id}{$new_linked_feat_id}) < $feat_hash{$new_linked_feat_id}->{'full'})) {
		    print STDERR "$new_linked_feat_id failed minimal linkage test\n" if ($DEBUG);
		    next; #link less than threshold
		}
		push (@linked_array, $new_linked_feat_id);
		#Remove links we have already used to add fragments
		print STDERR "frameshift weight: $con_comp_frames{$linked_feat_id}{$new_linked_feat_id}\n" if ($DEBUG);
		print STDERR "new_linked_feat_id = $new_linked_feat_id :full $feat_hash{$new_linked_feat_id}->{'full'}\n" if ($DEBUG);
		delete $con_comp_frames{$linked_feat_id}{$new_linked_feat_id};
		delete $con_comp_frames{$new_linked_feat_id}{$linked_feat_id};
	    }

	    #need to delete all references to this $linked_feat_id since it is a protein fragment we don't want to cluster
	    #need to rebuild %genome_hash and %TagByPointer after deleting all of the protein fragments
	    #do we want to update orf_counter to subtract out the frameshifts?
	    #$orf_counter{$FeatnameLookupTag_hash{$linked_feat_id}}->{'raw'}--;  # decrement the total orf counter when removing a fragment
	    #$orf_counter{$FeatnameLookupTag_hash{$linked_feat_id}}->{'used'}--;  # decrement the used orf counter when removing a fragment

	    #delete matches between the deprecated protein fragments and other proteins
	    foreach my $relate_id (keys %{ $relationship_hash{$linked_feat_id} }) {#remove symmetric relationships
		delete $relationship_hash{$relate_id}{$linked_feat_id};
	    }
	    delete $relationship_hash{$linked_feat_id};

	    my $qtag = $FeatnameLookupTag_hash{$linked_feat_id};
	    foreach my $stag (keys %{ $Qbytaghash{$linked_feat_id} } ) { #delete references symmetrically
		    foreach my $sid (keys %{ $Qbytaghash{$linked_feat_id}{$stag} })  {
			delete $Qbytaghash{$sid}{$qtag}{$linked_feat_id};
		    }
	    }
	    delete $Qbytaghash{$linked_feat_id};

	    delete $feat_hash{$linked_feat_id};
	    delete $Tagbyfeatnamehash{$qtag}{$linked_feat_id};
	    delete $FeatnameLookupTag_hash{$linked_feat_id};
	    delete $AssemblyLookup_hash{$linked_feat_id};
	    delete $TagByPointer{$linked_feat_id};

	}
	print OUTFRAMESHIFT "\n";
	    
	$prev_tag = $cur_tag;
	$prev_asmbl_id = $cur_asmbl_id;
    }
    close (OUTFRAMESHIFT);

    #rebuild %genome_hash and %TagByPointer after deleting all of the protein fragments
    &rebuild_genome_hash;
    return;
}
			    
sub make_tables {

############ Generate a match table of bidirectional best hits ###############
    
    open (OUTVENN, ">$outprefix" . "matchtable.txt") if ($vennfile);
    open (OUTZEROONE, ">$outprefix" . "matchtable_0_1.txt") if ($vennfile);
    open (OUTVENNID, ">$outprefix" . "matchtable_id.txt") if ($vennfile);
    open (OUTID, ">$outprefix" . "id.txt") if ($vennfile);
    open (MICRO, ">$outprefix" . "micro.txt") if ($microarray == 1);
    open (PAULSEN, ">$outprefix" . "hits.txt") if ($hitsfile == 1);
    open (READ, ">$outprefix" . "BSR.txt") if ($normalizefile == 1);
    foreach my $cluster_id (1 .. $#cluster_size)  { # output in previously computed order of the clusters
	my $cluster_member_first = ${ $cluster_for_number[$cluster_id] }[0];
	my $genome_tag_index = $cluster_member_first->{'tag'};
	my $query_featname = $cluster_member_first->{'id'};
	my $i = $genome_tag_index + 1; #used to know how many columns to skip for some tables
	my $subject_tag_index;
	my $subject_featname = "";
	for ($subject_tag_index = 0; $subject_tag_index < $genome_tag_index; $subject_tag_index++) {
	    &write_files("3", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname, $cluster_id);
	}
	# We are at $genome_tag now so we are looking at self comparisons, generate output and move on
	&write_files("1", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname, $cluster_id);#$subject_featname has not been defined yet here!!
	$subject_tag_index++; #go past $genome_tag
	# Now iterate through sorted (by subject_tag_index) cluster array to print out members of cluster
	foreach my $cluster_member ( @{ $cluster_for_number[$cluster_id] } ) {
	    my $cluster_member_tag_index = $cluster_member->{'tag'};
	    if ($cluster_member_tag_index < $subject_tag_index) {
		if ($cluster_member->{'id'} ne $query_featname) {
		    die ("ERROR!!! cluster for $query_featname should have been output sooner!\n");
		}
		next; #skip $genome_tag member we already output for
	    }
	    while ($cluster_member_tag_index > $subject_tag_index) {
		#output blanks for this genome
		&write_files("3", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname, $cluster_id);
		$subject_tag_index++; #go to next genome
	    }
	    $subject_featname = $cluster_member->{'id'};
	    &write_files("2", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname, $cluster_id);
	    $subject_tag_index++; #go to next genome
	}
	while ($subject_tag_index <= $#tag_array) {
	    &write_files("3", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname, $cluster_id);
	    $subject_tag_index++; #go to next genome
	}
	&write_files("4", $i, $subject_tag_index + 1, $query_featname, $tag_array[$subject_tag_index], $subject_featname, $cluster_id); #add newlines to files
    }
    close (OUTVENN) if ($vennfile);
    close (MICRO) if ($microarray == 1);
    close (PAULSEN) if ($hitsfile == 1);
    close (READ) if ($normalizefile == 1);
    close (OUTVENNID) if ($vennfile);
    close (OUTID) if ($vennfile);
}

sub read_clusters {

    open (my $clusterfile, "<", $cluster_input_file) || die ("ERROR: can't open file $cluster_input_file\n");
    my $count = 0;
    while (<$clusterfile>)  {
	if ( /^\s*$/ ) { #line has only whitespace so ignore
	    next;
	}
	chomp;
	$count++;
	my $cluster_num;
	my @cluster_members;
	($cluster_num, @cluster_members) = split(/\t/,$_);
	if ($count != $cluster_num) {
	    die ("ERROR: clusters are not consecutively numbered from 1 to N for file $cluster_input_file ($count:$cluster_num)\n");
	}
	$cluster_size[$cluster_num] = 0;
	my %seen_tag_already = (); # make sure only one cluster member per genome
	my $merged = [];
	my $prev_tag_index = -1; # this is less than any real tag index
	my $needs_sorting = 0; # if genome columns of the clusters file are in the same order as the genome identifier file we do not need to sort
	foreach my $member (@cluster_members)  { # add all cluster members to the cluster
	    if ($member eq '') {
		next; # skip null columns
	    }
	    my $tag;
	    my $tag_index;
	    my ($qid) = split(/\s+/, $member);  # split $member on space in case there is anything like (%identitiy) after the qid
	    if ($qid eq '') {
		next; # skip columns that have nothing but spaces
	    }
	    if (!defined ($tag = $FeatnameLookupTag_hash{$qid})) {
		die("ERROR!!! for cluster $cluster_num cluster member $qid was not associated with a genome indentifier in the input attribute file\n");
	    }
	    if (defined $seen_tag_already{$tag}) {
		die("ERROR!!! for cluster $cluster_num the genome identifier $tag for cluster member $qid was the same as cluster member $seen_tag_already{$tag} but only one cluster member per genome is allowed\n");
	    } else {
		$seen_tag_already{$tag} = $qid;
	    }
	    if (!defined ($tag_index = $TagIndex{$tag})) {
		die("ERROR!!! genome identifier $tag for cluster member $qid was not associated with an index into the genome identifier array\n");
	    }
	    if ($tag_index <= $prev_tag_index) {
		$needs_sorting = 1; # now we need to sort
	    }
	    $prev_tag_index = $tag_index;
	    push(@{ $merged }, { 'id' => $qid, 'tag' => $tag_index });
	    $clusters{$qid} = $merged;
	    $cluster_number{$qid} = $cluster_num;
	    $cluster_size[$cluster_num]++;
	}
	$cluster_for_number[$cluster_num] = $merged;
	if ($needs_sorting) {
	    @{ $merged } = sort {$a->{'tag'} <=> $b->{'tag'}} ( @{ $merged } );
	}
    }
    close($clusterfile);

    my $changed = 0;
    foreach my $qid (keys %feat_hash)  { # go through all featnames
	if (!defined $cluster_number{$qid})  {# for featnames in the attribute file but not in clusters (probably frameshift fragments) cleanup hashes to ignore these proteins
	    $changed++;
	    delete $feat_hash{$qid};
	    delete $FeatnameLookupTag_hash{$qid};
	    delete $AssemblyLookup_hash{$qid};
	    delete $TagByPointer{$qid};
	}
    }
    if ($changed) {
	print STDERR "WARNING: $changed features in the attribute file were not in the clusters file - probably frameshifted fragments\n";
	&rebuild_genome_hash; #need to remove the deleted features/genes/proteins
    }

    return();

}

sub print_report  {

    my $key = "";
    open (OUT, ">$outprefix" . "report.txt");
    print OUT " e-value cut-off: <= $evalue\n";
    print OUT "percent identity: >= $percentid\n";
    print OUT " length of match: >= $min_hit_length\n";
    print OUT "input .btab file: $btabpath/$btabfile\n\n";
    foreach $key (@tag_array) {
	print OUT "Raw $key = $orf_counter{$key}->{'raw'}\n";
	print OUT "Used $key = $orf_counter{$key}->{'used'}\n";
    }
    close (OUT);
    return();

}

sub option_help {

   system("clear");
   print STDERR <<_EOB_;
$prog:

           Pan-genome Ortholog Clustering Tool, a heuristic computer program for 
           pan-genomic analysis of closely related prokaryotic species or strains.

Copyright (C) 2011-2015  The J. Craig Venter Institute (JCVI).  All rights reserved

License:   This program is free software: you can redistribute it and/or modify
           it under the terms of the GNU General Public License as published by
           the Free Software Foundation, either version 3 of the License, or
           (at your option) any later version.

           This program is distributed in the hope that it will be useful,
           but WITHOUT ANY WARRANTY; without even the implied warranty of
           MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
           GNU General Public License for more details.

           You should have received a copy of the GNU General Public License
           along with this program.  If not, see <http://www.gnu.org/licenses/>.

Citation:  Derrick E. Fouts, Lauren Brinkac, Erin Beck, Jason Inman, and Granger Sutton 
           (2012) "PanOCT: Automated Clustering of Orthologs using Conserved Gene 
           Neighborhood for Pan-Genomic Analysis of Bacterial Strains and Closely Related 
           Species" Nucleic Acids Res. 40(22):e172.

  Usage: $prog <options>
Example: $prog -t example_blast.txt -f example_tags.txt -g example.gene_att -P example.pep -S Y -L 1 -M Y -H Y -V Y -N Y -F 1.33 -G y -c 0,25,50,75,100 -T
         or
         $prog -t example_blast.txt -f example_tags.txt -g example.gene_att -P example.pep -i 20 -F 1.33 -N Y -M Y -H Y
         or
	 if clusters have already been calculated (from another tool or the hierarchical clustering version of the pan-genome pipeline available for download in the svn of PanOCT) and you wish to generate the PanOCT output files:
	 $prog -R example_clusters.txt -f example_tags.txt -g example.gene_att -P example.pep
Version: $version
Options:
     -h: print this help page
     -c: argument is a comma separated list of numbers between 0-100 (0,50,75,100 might be a good choice) each number represents the percent of genomes
         needed for a cluster to be considered a core cluster. For each number two files are generated: one for core clusters and one for noncore
         clusters. The core file shows the order and orientation of core clusters with respect to each other. The noncore file shows the order and
         orientation of noncore clusters with respect to core clusters.
     -e: argument is a comma separated list of numbers between 0-100 each number represents the percent of genomes needed for a cluster to be considered
         a core cluster. For each number BSR distance matrices are computed using only the "core" clusters. [DEFAULT = 0,100]
     -d: no argument, generates a mutlifasta file of cluster centroids [on by DEFAULT]
     -T: no argument, prints cluster numbers as the first column in all table files [off by DEFAULT]
     -W: window size on either side of match to use CGN [DEFAULT = 5] must be between [1,20]
     -b: base directory path [DEFAULT = PWD]
     -p: path to btab file [DEFAULT = base directory]
     -t: name of btab (wublast-style or ncbi -m8 or -m9) input file [REQUIRED]
     -f: file containing unique genome identifier tags [REQUIRED]
     -g: gene attribute file (asmbl_id<tab>protein_identifier<tab>end5<tab>end3<tab>annotation<tab>genome_tag)
     -P: name of concatinated .pep file [REQUIRED to calc protein lengths]
     -Q: path to .pep file [DEFAULT = base directory]
     -i: aa % identity cut-off [DEFAULT = 35.0]
     -I: aa % identity cut-off for frameshift detection [DEFAULT = 35.0]
     -E: E-value [DEFAULT = 0.00001]
     -L: Minimum % match length [DEFAULT = 1]
     -H  Want to create a file with the id and annotation of each protein in a cluster?  (y)es or (n)o [DEFAULT = NO]
     -V: Want to create an ortholog matchtable?  (y)es or (n)o [DEFAULT = YES]
     -N: Want to create a normalized BLAST score file?  (y)es or (n)o [DEFAULT = NO]
     -M: Want to create microarray-like data for normalized BLAST scores?  (y)es or (n)o [DEFAULT = NO]
     -H: Want to create a table of hits (y)es or (n)o [DEFAULT = NO]
     -G: Want to create a normalized BLAST score histogram file containing a row for each genome and two rows for each pair of genomes?  (y)es or (n)o [DEFAULT = NO]
     -C: Want to create a file with the number of matches within clusters and statistics on protein length per cluster?  (y)es or (n)o [DEFAULT = YES]
     -A: Want to create a file grouping ortholog clusters which appear to be close paralogs?  (y)es or (n)o [DEFAULT = YES]
     -U: Want to create a file with ids of proteins which appear to be fragments or fusions based on clustering?  (y)es or (n)o [DEFAULT = YES]
     -B: Want to create a file with pairwise similarity scores used for clustering?  (y)es or (n)o [DEFAULT = YES]
     -S: Want to use strict criteria for ortholog determination?  (y)es, (m)intermediate or (n)o [DEFAULT = YES]
     -F: Deprecate shorter protein fragments when protein is split due to frameshift or other reason
         Takes an argument between 1.0 and 2.0 as a length ratio test - recommended value is 1.33 [DEFAULT = off]
     -a: Number of amino acids at the beginning or end of a match that can be missing and still be
         considerd a full length match - must be between 0 and 100 - [DEFAULT = 20]
     -s: Number of blast matches needed to confirm a protein fragment/frameshift [DEFAULT = 1]
     -R: argument is a file name for clusters generated by PanOCT or some othr program in PanOCT format, PanOCT will not recluster but simply generate some of its output files
     -D: DEBUG MODE (DEFAULT = off)

 Authors: Derrick E. Fouts, Ph.D. and Granger Sutton, Ph.D.
 Date: December 21, 2004; last revised August 18, 2015
 Note: For detailed descriptions of input and output files, please read the README.txt file found in this distribution

_EOB_
    exit;
}

########################################  M A I N  #####################################################
open (my $parameterfile, ">$outprefix" . "parameters.txt");
&print_parameters ($parameterfile);
close ($parameterfile);
print STDERR "fetching tags to search\n";
&get_tags;
print STDERR "Gathering protein sequence information from $pep_file\n";
&get_protein_info;
print STDERR "gathering gene attributes\n";
&get_gene_att;
if ($process_clusters) {
    &read_clusters;
    if ($compute_adjacency) {
        print STDERR "Calculating cluster adjacency ...\n";
        &calc_adjacency;
    }
    print STDERR "Calculating pairwise cluster distances ...\n";
    &calc_pairwise_cluster_distances;
} else {
    print STDERR "getting data from .btab file\n";
    &select_max_scores_from_btab;
    &select_data_from_btab;
    print STDERR "Calculating the Blast Score Ratio (BSR) ...\n";
    &calc_BSR;
    if ($frameshiftfile) {
        print STDERR "Determining frameshifts / protein fragments ...\n";
        &frameshifts;
    }
    print STDERR "Calculating reciprocal best hits (RBHs) ...\n";
    &calc_bibest;
    print STDERR "Calculating cliques of RBHs ...\n";
    &calc_cliques;
    print STDERR "Calculating anchors ...\n";
    &calc_synteny ($ANCHOR_window_size);
    &calc_anchors;
    print STDERR "Extending anchors ...\n";
    &extend_anchors;
    print STDERR "Calculating conserved gene neighborhood (CGN) scores ...\n";
    &calc_synteny ($CGN_window_size);
    print STDERR "Calculating recipricol best CGN scores ...\n";
    &calc_synbibest;
    if ($histogramfile || $strict_orthos) {
        print STDERR "Calculating histograms ...\n";
        &calc_score_histograms;
    }
    print STDERR "Calculating clusters ...\n";
    &calc_clusters;
    &calc_cluster_numbers;
    if ($compute_centroids) {
        print STDERR "Calculating cluster centroids ...\n";
        &calc_cluster_centroids;
    }
    print STDERR "Calculating cluster weights ...\n";
    &calc_cluster_weights; # These are numbers of good matches between clusters and within a cluster - also outputs paralog clusters
    if ($compute_adjacency) {
        print STDERR "Calculating cluster adjacency ...\n";
        &calc_adjacency;
    }
    if ($pairwisefile)  {
        print STDERR "Printing pairwise matches within clusters ...\n";
        &print_pairwise_cluster_matches;
    }
    print STDERR "Calculating BSR distances ...\n";
    &calc_all_BSR_distances;
    print STDERR "Calculating pairwise cluster distances ...\n";
    &calc_pairwise_cluster_distances;
    print STDERR "building table for microarray software\n" if ($microarray);
    print STDERR "building table of hits\n" if ($hitsfile);
    print STDERR "generating match table!\n" if ($vennfile);
    &make_tables;
    &print_report;
}
exit(0);
