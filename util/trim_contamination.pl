#!/usr/local/bin/perl -w

#################################################################################
#										#				
# Name	      : trim_contamination.pl						#	
# Version     : 1.0								#
# Project     : GenBank Submissions						#
# Description : Script to trim contanimated sequences out of FASTA and table files	#
# Author      : Shaun Adkins							#
# Date        : April 8, 2013							#
#										#
#################################################################################

use strict;
# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
# Prints usage message from embedded pod documentation
use Pod::Usage;
# Include packages here 
use File::Basename;
use FindBin qw($Bin);
use List::Util qw(first);

###########
# GLOBALS #
###########
my %cmdLineArgs = ();
# Log file handle;
my $logfh;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);
my $outdir = {};
my $output_file = {};
my $genome_ref = {};
my $contaminated_contigs = {};
my @cols = ();
my ($line, $v, $file, $head, $m, $contig);
my $flag = 0;
my $seq = "";
my $table = 0;

################
# MAIN PROGRAM #
################
GetOptions(\%cmdLineArgs,
	   'contaminants|c=s',
	   'input_file|f=s',
	   'output_dir|o=s',
	   'table|t',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($cmdLineArgs{'help'});

checkCmdLineArgs();

$file = $cmdLineArgs{'input_file'};
my @file_split = split("/", $file);	# Just the basename
$outdir = $cmdLineArgs{'output_dir'};
system("mkdir $outdir") if (! -e $outdir);	
$output_file = $outdir . "/" . "new_" . $file_split[-1];	# Create our output file

$contaminated_contigs = parseContaminants($cmdLineArgs{'contaminants'});

if ($table) {
	process_table($file, $output_file);
} else {
	process_fasta($file, $output_file);
}

###############
#       SUBROUTINES      #
###############

sub process_fasta {
	my ($file, $fsa_file) = @_;
	$genome_ref = {};
	$genome_ref = readSeq($file);
	open(FO, "> $fsa_file") or printLogMsg($ERROR, "Could not open file $fsa_file for writing. Reason: $!");
	foreach $contig (keys %$genome_ref) {
		if(exists($contaminated_contigs->{$contig})) {
			# Start and end coordinates assume sequence start at position 1
			my ($start, $end) = split(/\.\./, $contaminated_contigs->{$contig});
			my $new_seq = $genome_ref->{$contig};
			$new_seq = substr($new_seq, 0, $start-1) . substr($new_seq, $end);
			#print "New Seq for $contig is from 0.." . ($start-1) . " and from $end to the end\n";
			print FO ">$contig\n$new_seq\n";
		} else {
			print FO ">$contig\n$genome_ref->{$contig}\n";
		}
	}
	close(FO);
	#system("$Bin/cleanFasta.pl $fsa_file");	

}

# Checking table gene coordinates against contamination coordinates
# Make adjustments where necessary.
# If contamination affects coordinates in the middle or to the left of the gene, then the gene and all following genes will have their coordinates shifted to the left to close the created gap.
# For example if the gene is from 1..80 and the contamination is from 20..30, then the gene will be split 1..>19 and >20..70. 
sub process_table {
    my ($file, $tbl_file) = @_;
    my $tbl_hash = readTbl($file);
    
    foreach my $contig (keys %{$contaminated_contigs}) {
    	my ($start, $end) = split(/\.\./, $contaminated_contigs->{$contig});
    	my $shift_flag = 0;
    	print"$contig\n";
    	
    	# sort gene arrays by coordinates via Schwartzian transformation
    	if (%{ $tbl_hash->{$contig}->[0] }) {
    		@{$tbl_hash->{$contig}} = map {$_->[0]}
    							   	   sort {$a->[1] <=> $b->[1]} 
    							   	   map {$_->{'left'} =~ /[<>]?(\d+)/; [$_, $1] }
    							   	   @{$tbl_hash->{$contig}};
    	}
    	
    	foreach my $gene (@{$tbl_hash->{$contig}}) {
    	    next if ! (%{$gene});	# In rare cases where there were no entries for a given feature;

    	    my ($first, $second) = ($gene->{'left'}, $gene->{'right'});
    	    $first =~ s/[<>]//;	#get rid of partial signs
    	    $second =~ s/[<>]//;
	    $gene->{'shift'} = $shift_flag;	# How many positions should we shift this gene to the left.
	    my $pid = $gene->{'protein_id'};	# How many positions should we shift this gene to the left.
	    
    	    # Determine where contamination overlaps and adjust coordinates
    	    if ($start > $first && $end < $second) { 
    	        # The contamination splits the gene
    	        my $position = first { ${$tbl_hash->{$contig}}[$_]->{'locus_tag'} eq $gene->{'locus_tag'} } 0..$#{$tbl_hash->{$contig}}; 	# Get position of gene in array
     	        my %new_split = %{$gene};
     	        if (defined $new_split{'protein_id'}) {
     	            $gene->{'protein_id'} = $gene->{'protein_id'} . "_A";
     	            $new_split{'protein_id'} = $new_split{'protein_id'} . "_B";
     	        }
     	        $gene->{'right'} = ($gene->{'revcom'}) ? "<" . ($start-1) : ">" . ($start -1);
     	        $new_split{'left'} = ($gene->{'revcom'}) ? ">" . ($end+1) : "<" . ($end+1);
     	        $new_split{'right'} = $second;
     	        $shift_flag += ($end - $start + 1 );	# If a previous shift occurred in this gene, accumulate shiftting additions
     	        splice @{$tbl_hash->{$contig}}, $position+1, 0, \%new_split;
     	        $second = ($start -1);	# Since we are still processing our "A" gene, need to adjust $second 
     	                    	            
    	    } else {
    	        # The contamination overlaps one of the boundaries of the gene
    	        if ($end >= $first && $end < $second) { # The overlap extends beyond the left
    	            $first = $start;
     	            $second = $second - ($end - $start + 1);   	                
    	            $gene->{'left'} = ($gene->{'revcom'}) ? ">" . $first : "<" . $first;
     	            $shift_flag += ($end - $start + 1);
    	        }
    	        elsif ($start > $first && $start <= $second) {  # The overlap extends beyond the right
    	            $second = $start -1;
    	            $gene->{'right'} = ($gene->{'revcom'}) ?  "<" . $second : ">" . $second;
    	            $shift_flag += ($end - $start + 1);
    	        }
    	        elsif ($first > $start && $first > $end) {  # The trimming happens entirely on the left side
    	            $first = ($first - ($end - $start + 1));
    	            $second = ($second - ($end - $start + 1));
    	        }
    	    }
    	    
    	    # applying shift in coordinate position.
    	    $first = $first - $gene->{'shift'};
    	    $second = $second - $gene->{'shift'};
    	    $gene->{'left'} =~ s/\d+/$first/;
    	    $gene->{'right'} =~ s/\d+/$second/;
    	    # setting codon_start position.... if partial start is present... delete codon_start key
    	    if ($gene->{'left'} =~ /</) {
    	    	delete $gene->{'codon_start'};
    	    } else {
    	    	$gene->{'codon_start'} = write_codon_pos($first);
    	    }
    	}
    }
    
    write_tbl ($tbl_file, $tbl_hash);
}

sub readTbl {
    my $file = shift;
    
    my %tbl_hash;
    my @gene;
    my %new_gene;
    my $contig;
    open (FR, "<$file") or printLogMsg($ERROR, "Could not open file $file for reading. Reason:  $!");
    
    while ($line = <FR>) {
        chomp($line);
        next if ($line =~ /^\s*$/);
        if($line =~ /^>Feature\s+(\S+)/) {
            # Add last genes and contigs before processing new one
            if (defined $contig) {
                my %passed_gene = %new_gene;
                push @gene, \%passed_gene;
                my @g = @gene;
            	$tbl_hash{$contig} = \@g;            
            }
            # Process new contig
	    $contig = $1;
	    @gene = ();	# Reset array for new sets of genes
	    %new_gene = ();	# Reset hash to store gene info
        } elsif ($line =~ /^([<>]?\d+)\s+([<>]?\d+)\s+gene/) {
            if (%new_gene) {
                my %passed_gene = %new_gene;
		push @gene, \%passed_gene;	#if contig has multiple genes, start storing into array
	    }
            %new_gene = ();
            $new_gene{'left'} = $1;
            $new_gene{'right'} = $2;
            my ($first, $second) = ($new_gene{'left'}, $new_gene{'right'});
    	    $first =~ s/[<>]//;	#get rid of partial signs
    	    $second =~ s/[<>]//;
            my $revcom = is_revcom($first, $second);
            if ($revcom) {
                my $temp = $new_gene{'left'};
                $new_gene{'left'} = $new_gene{'right'};
                $new_gene{'right'} = $temp;
                $new_gene{'revcom'} = 1;
            }
            
        } elsif ($line =~/^([<>]?\d+\s+){2}(CDS|tmRNA|tRNA|rRNA|mRNA)/) {
            $new_gene{'type'} = $2;
        } else {
            #attempting to capture every other line
            if ($line =~/^\s+([a-zA-Z_]+)\s+(.+)/ ){
        	my $part = $1;
        	$new_gene{$part} = $2;
            } elsif ($line =~ /^\s+pseudo/){
            	$new_gene{'pseudo'} = 1;
            } else {
        	print "Line $line from $contig is a bit wonky\n";
            }
        }
    }
    
     push @gene, \%new_gene;
     $tbl_hash{$contig} = \@gene;     
     close FR;
     
     return \%tbl_hash;
}

sub readSeq {
	my ($file) = @_;
	my ($head, $seq, $line);
	my %genome = ();
	open(FR, "< $file") or printLogMsg($ERROR, "Could not open file $file for reading. Reason: $!");
	while($line = <FR>) {
		chomp($line);
		next if ($line =~ /^\s*$/);
		if($line =~ /^>(.+)/) {
			if($seq) {
				$genome{$head} = $seq;	
			}
			$head = $1;
			$head =~ s/\s+$//;	
			$seq = "";	
		} else {
			$seq .= $line;
		}
	}
	$genome{$head} = $seq;
	close(FR);
	return(\%genome);
}

sub write_tbl {
	my ($tbl_file, $tbl_hash) = @_;
	
	open TBL, ">$tbl_file" or printLogMsg($ERROR, "Cannot open $tbl_file for writing. Reason $!");
	open(CH, ">$outdir/contamination.tbl") or printLogMsg($ERROR, "Could not open file $outdir/contamination.tbl for writing. Reason: $!");

	#Schwartzian transform to sort by contig number
	my @sorted = map {$_ -> [0]}							#Return original key value
				   sort { $a->[1] <=> $b->[1] || $a cmp $b}	#Sort by contig number
	                           map {$_ =~ /\.(\d+)$/; [$_, $1]}			# Map regex to anon arrayref
	                           keys %$tbl_hash;
	
	# Now to write the table
	foreach my $contig (@sorted) {
	 	print TBL ">Feature $contig\n";
	 	print CH ">Feature $contig\n";
	 	my $contam_flag = 0;
	 	$contam_flag = 1 if (exists $contaminated_contigs->{$contig});
	 	
    		foreach my $gene (@{$tbl_hash->{$contig}}) {
    		    	next if ! (%{$gene});	# In rare cases where there were no entries for a given feature;
	 		print TBL $gene->{'left'}, "\t", $gene->{'right'}, "\tgene\n" if (! defined $gene->{'revcom'});
	 		print TBL $gene->{'right'}, "\t", $gene->{'left'}, "\tgene\n" if (defined $gene->{'revcom'});	 	    	
	 		print TBL "\t\t\tlocus_tag\t", $gene->{'locus_tag'}, "\n";
	 		print TBL "\t\t\tgene\t", $gene->{'gene'}, "\n" if (defined $gene->{'gene'});
	 		print TBL $gene->{'left'}, "\t", $gene->{'right'}, "\t", $gene->{'type'}, "\n" if (! defined $gene->{'revcom'});
	 		print TBL $gene->{'right'}, "\t", $gene->{'left'}, "\t", $gene->{'type'}, "\n" if (defined $gene->{'revcom'});
	 		
	 		if ($contam_flag){
	 			print CH $gene->{'left'}, "\t", $gene->{'right'}, "\tgene\n" if (! defined $gene->{'revcom'});
	 			print CH $gene->{'right'}, "\t", $gene->{'left'}, "\tgene\n" if (defined $gene->{'revcom'});	 	    	
	 			print CH "\t\t\tlocus_tag\t", $gene->{'locus_tag'}, "\n";
	 			print CH "\t\t\tgene\t", $gene->{'gene'}, "\n" if (defined $gene->{'gene'});
	 			print CH $gene->{'left'}, "\t", $gene->{'right'}, "\t", $gene->{'type'}, "\n" if (! defined $gene->{'revcom'});
	 			print CH $gene->{'right'}, "\t", $gene->{'left'}, "\t", $gene->{'type'}, "\n" if (defined $gene->{'revcom'});
	 		}
	 		
	 	 	delete $gene->{'left'}; delete $gene->{'right'};
	 	 	delete $gene->{'type'}; delete $gene->{'gene'}; delete $gene->{'locus_tag'};
	 	 	
	 	 	foreach my $feat (keys %{$gene}) {
	 	 		next if ($feat eq "revcom");	# skip hash keys that was created to aid in processing
	 	 		next if ($feat eq "shift");	
	 	 		if ($feat eq "pseudo") {	# Pseudomolecules had no values
	 	 			print TBL "\t\t\tpseudo\n";
	 	 			print CH "\t\t\tpseudo\n" if ($contam_flag);
	 	 		} else {
	 	 			print TBL "\t\t\t", $feat, "\t", $gene->{$feat}, "\n";
	 	 			print CH "\t\t\t", $feat, "\t", $gene->{$feat}, "\n" if ($contam_flag);
	 	 		}
	 	 	}
	 	}
	}
	
	close TBL;
	close CH;
}

####################################################################################################################################################

sub parseContaminants {
	my ($filename) = @_;
	my $line;
	my %contaminants = ();
	open(FR, "< $filename") or printLogMsg($ERROR, "Could not open file $filename for reading. Reason: $!");
	while($line = <FR>) {
		chomp($line);
		next if ($line =~ /^\s*$/);
		# contig_name	seq_length	start..stop	adaptor
		my @arr = split(/\s+/, $line);
		$contaminants{$arr[0]} = $arr[2];
	}
	close(FR);
#	foreach my $t (keys %contaminants) {
#		print "$t:\t@{$contaminants{$t}}\n";
#	}
	return(\%contaminants)
}

#############################################################################

sub write_codon_pos {
	my $start = shift;
	my $num = $start%3;
	$num = 3 if ($num == 0);
	return $num;
}

####################################################################################################################################################

sub is_revcom {
	my ($first, $second) = @_;
	return 1 if ($second < $first);
	return 0 if ($first > $second);
}

#############################################################################

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
	my @required = qw(contaminants input_file output_dir);
	foreach my $option(@required) {
		if(!defined($cmdLineArgs{$option})) {
			 printLogMsg($ERROR,"ERROR! : Required option $option not passed");
		}
	}
	if(exists($cmdLineArgs{'log'})) {
		open($logfh, "> $cmdLineArgs{'log'}") or die "Could not open $cmdLineArgs{'log'} file for writing.Reason : $!\n"
	}
	
	$table = 1 if ($cmdLineArgs{'table'});
}

####################################################################################################################################################

# Description   : Used to handle logging of messages(errors and warnings) during the execution of the script
# Parameters    : level = can be ERROR, WARNING or INFO
#		  msg   = msg to be printed in the log file or to STDERR
# Returns       : NA
# Modifications : 

sub printLogMsg {
	my ($level, $msg) = @_;
	if( $level <= $DEBUG ) {
		print STDERR "$msg\n";
		print $logfh "$msg\n" if(defined($logfh));
		die "\n" if($level == $ERROR);
	}	
}

__END__

#####################
# POD DOCUMENTATION #
#####################

=head1 NAME

# Name of the script and a 1 line desc

=head1 SYNOPSIS

# USAGE : 

	parameters in [] are optional

=head1 OPTIONS



=head1 DESCRIPTION



=head1 INPUT



=head1 OUTPUT



=head1 AUTHOR

	Shaun Adkins
	Bioinformatics Software Engineer I
	Institute for Genome Sciences
	School of Medicine
	University of Maryland
	Baltimore, MD - 21201
	sadkins@som.umaryland.edu

==cut
