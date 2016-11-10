#!/usr/local/bin/perl -w

#########################################################################################
#											#
# Name	      : gene_symbol_mod_tbl.pl							#
# Version     : 1.0									#
# Project     :	GenBank Submisison Pipeline						#
# Description : Script to remove upper case, duplicate gene symbols as well replace 	#
#		existing gene symbols in tbl with new ones provided in the rules file	#
# Author      : Sonia Agrawal								#
# Date        : February 04, 2014							#
#											#
#########################################################################################

use strict;

# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

# Prints usage message from embedded pod documentation
use Pod::Usage;

# Include packages here
use DBI;
use File::Basename;

#############
# CONSTANTS #
#############

###########
# GLOBALS #
###########
my %hCmdLineArgs = ();

# Log file handle;
my $logfh;
my ( $ERROR, $WARN, $DEBUG ) = ( 1, 2, 3 );
my %hFTP  = ();
my @aMeta = ();

################
# MAIN PROGRAM #
################
GetOptions(
    \%hCmdLineArgs, 'input_file|i=s', 'output_dir|o=s', 'tbl_file|t=s', 'log|l=s',
    'help|h'
) or pod2usage();

pod2usage( { -exitval => 0, -verbose => 2, -output => \*STDERR } )
  if ( $hCmdLineArgs{'help'} );

checkCmdLineArgs( \%hCmdLineArgs );

readInput( $hCmdLineArgs{'input_file'}, \@aMeta );

if (   defined( $hCmdLineArgs{'tbl_file'} )
    && ( -e $hCmdLineArgs{'tbl_file'} )
    && ( -s $hCmdLineArgs{'tbl_file'} )
    && defined( $hCmdLineArgs{'output_dir'} ) )
{
    fixGeneSymbolsInTbl( $hCmdLineArgs{'tbl_file'},
        $aMeta[3], $hCmdLineArgs{'output_dir'} );
} else {
    printLogMsg( $ERROR,
        "ERROR : A valid .tbl file is needed to fix the gene symbols in .tbl file."
    );
}

###############
# SUBROUTINES #
###############

####################################################################################################################################################
# Description   : Used to remove duplicate and beginning with upper-case gene symbols. Also replaces gene symbols based on rules file if provided.
# Parameters    : sTblFile, sRFile, sOutDir
#		  sTblFile - path to original table file to be fixed
#		  sRFile - path to rules file with gene symbols to be replaced
#		  sOutDir - path to output directory to create the fixed .tbl file
# Returns       : NA
# Modifications :

sub fixGeneSymbolsInTbl {
    my ( $sTblFile, $sRFile, $sOutDir ) = @_;
    my $sSubName   = ( caller(0) )[3];
    my %hRules     = ();
    my %hExistSyms = ();
    my ( $fhRead, $fhWrite );
    my ( $sLine, $sSymb, $sFBase, $sFDir, $sFExt, $sCorrectedTbl, $sGeneSymb );
    my $nCnt = 0;

    if ( defined($sRFile) && ( -e $sRFile ) && ( -s $sRFile ) ) {
        readRulesFile( $sRFile, \%hRules );
    } else {
        printLogMsg( $WARN,
            "WARNING : $sSubName :: Gene symbols rules does not exist or is empty. Thus only duplicate and upper case gene symbols  will be removed from $sTblFile file"
        );
    }

    open( $fhRead, "< $sTblFile" )
      or printLogMsg( $ERROR,
        "ERROR : $sSubName :: Could not open file $sTblFile for reading. Reason : $!"
      );
    while ( $sLine = <$fhRead> ) {
        chomp($sLine);
        next if ( $sLine =~ /^\s+$/ );
        if ( $sLine =~ /\d+\s+\d+\s+gene/ ) {
            if ($sSymb) {
                if ( exists( $hExistSyms{$sSymb} ) ) {
                    $hExistSyms{$sSymb} += 1;
                } else {
                    $hExistSyms{$sSymb} = 1;
                }
            }
            undef $sSymb;
        } elsif ( $sLine =~ /^\s+gene\s+(\S+)/ ) {
            $sSymb = $1;
        }
    }
    close($fhRead);

   # Create a list of duplicate and upper case gene symbols existing in tbl file
    foreach $sGeneSymb ( keys %hExistSyms ) {
        if ( ( $hExistSyms{$sGeneSymb} > 1 ) || ( $sGeneSymb =~ /^[A-Z]/ ) ) {
            $hRules{$sGeneSymb} = "";
        }
    }

    open( $fhRead, "< $sTblFile" )
      or printLogMsg( $ERROR,
        "ERROR : $sSubName :: Could not open file $sTblFile for reading. Reason : $!"
      );

    ( $sFBase, $sFDir, $sFExt ) = fileparse( $sTblFile, qr/\.[^.]*/ );

    if ( ( -d $sOutDir ) && ( -e $sOutDir ) ) {
        $sCorrectedTbl = $sOutDir . "/" . $sFBase . "_gs_corrected" . $sFExt;
        open( $fhWrite, "> $sCorrectedTbl" )
          or printLogMsg( $ERROR,
            "$sSubName :: Could not open file $sCorrectedTbl for reading. Reason : $!"
          );
    } else {
        printLogMsg( $ERROR,
            "ERROR : $sSubName :: Invalid output directory $sOutDir specified. Directory should exist"
        );
    }

    while ( $sLine = <$fhRead> ) {
        chomp($sLine);
        if ( $sLine =~ /^(\s+)gene(\s+)(\S+)/ ) {
            if ( exists( $hRules{$3} ) ) {
                $nCnt++;
                if ( $hRules{$3} eq "" || length( $hRules{$3} ) == 0 ) {
                    printLogMsg( $DEBUG,
                        "INFO : $sSubName :: Removing invalid gene symbol $3 from $sTblFile"
                    );
                    next;
                } else {
                    $sLine = $1 . "gene" . $2 . $hRules{$3};
                }
            }
        }
        if ( $sLine =~
            m/conserved hypothetical protein|conserved domain protein|conserved protein/
          )
        {
            printLogMsg( $DEBUG,
                "INFO : $sSubName :: Replacing $sLine with hypothetical protein in $sTblFile"
            );
            $sLine =~
              s/conserved hypothetical protein|conserved domain protein|conserved protein/hypothetical protein/gi;
        }
        if ( $sLine =~ m/tRNA-Pseudo/ ) {
            printLogMsg( $WARN,
                "WARNING : $sSubName :: tRNA-Pseudo issue encountered on line : $sLine"
            );
        }
        print $fhWrite $sLine . "\n";
    }
    printLogMsg( $DEBUG,
        "INFO : $sSubName :: Total number of gene symbols changed or deleted in tbl file $sTblFile is : $nCnt"
    );
    close($fhRead);
    close($fhWrite);
}

####################################################################################################################################################
# Description   : Used to read gene symbols rules file containing new values and old values to be replaced.
# Parameters    : sFile, phRules
#		  sFile - path to gene symbols rules file
#		  phRules - pointer to hash of rules to be populated
# Returns       : NA
# Modifications :

sub readRulesFile {
    my ( $sFile, $phRules ) = @_;
    my $sSubName = ( caller(0) )[3];
    my $fhRead;
    my ( $sLine, $sNewVal, $sOldVal );

    open( $fhRead, "< $sFile" )
      or printLogMsg( $ERROR,
        "ERROR : $sSubName :: Could not open file $sFile for reading. Reason : $!"
      );
    while ( $sLine = <$fhRead> ) {
        chomp($sLine);
        next if ( $sLine =~ /#/ );
        next if ( $sLine =~ /^\s+$/ );
        ( $sNewVal, $sOldVal ) = split( /\t|\s+/, $sLine, 2 );
        $phRules->{$sOldVal} = defined($sNewVal) ? $sNewVal : '';
    }
    close($fhRead);
}

####################################################################################################################################################
# Description   : Used to read input meta data file
# Parameters    : sFile, paMeta
#		  sFile - Path to input meta data file containing the information as described below
#		  paMeta - pointer to array to hold meta data
# Returns       : NA
# Assumptions	: The input meta data file sFile contains only a single row with information about a single database. This script is designed to be
#		  run interatively within an ergatis component. So it processes one DB at a time. The ergatis component itself could process
#		  several databases.
# Modifications :

sub readInput {
    my ( $sFile, $paMeta ) = @_;
    my $nI;
    my $sLine;
    my $fhRead;
    my $sSubName = ( caller(0) )[3];

    open( $fhRead, "< $sFile" )
      or printLogMsg( $ERROR,
        "ERROR : $sSubName :: Could not open $sFile file for reading. Reason : $!"
      );
    while ( $sLine = <$fhRead> ) {
        chomp($sLine);
        next if ( $sLine =~ /^#/ );
        next if ( $sLine =~ /^\s+$/ );
        @{$paMeta} = split( /\t/, $sLine );

# @meta : This script needs columns 0 and 3
#		[0] = Db name
#		[1] = NCBI locus tag
#		[2] = Path to rules file
#		[3] = Path to gene symbols file
#		[4] = Bioproject Id
#		[5] = Organism name
#		[6] = strain name
#		[7] = Organism serotype
#		[8] = Host
#		[9] = Date
#		[10]= Country
#		[11]= Assembly method
#		[12]= Genome Coverage
#		[13]= Sequencing platform
#		[14]= Contact person's name. Last name\sFirst name
#		[15]= Contact person's email address
#		[16]= Comma-separated author list. Each name should be Last name\sFirst name\sMiddle Initial
#		[17]= Title of the publication

        if ( length( $paMeta->[0] ) == 0 ) {
            printLogMsg( $WARN,
                "WARNING : $sSubName :: Database name missing for $sLine. Seems like you mean to fix the .tbl file only"
            );
        }
        if (   length( $paMeta->[3] ) == 0
            || !( -e $paMeta->[3] )
            || !( -s $paMeta->[3] ) )
        {
            printLogMsg( $WARN,
                "WARNING : $sSubName :: Missing gene symbols rules file in input file $sFile in column 4. Thus only duplicate and upper case gene symbols will be removed from either provided .tbl file or database."
            );
        }

    }
    close($fhRead);
}

####################################################################################################################################################
# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing.
# Parameters    : phCmdLineArgs
#		  phCmdLineArgs - reference to hash of command line arguments passed to the perl script
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
    my ($phCmdLineArgs) = @_;
    my $sOption;
    my @aRequired = ();
    my $sSubName  = ( caller(0) )[3];

    if ( exists( $phCmdLineArgs->{'log'} ) ) {
        open( $logfh, "> $phCmdLineArgs->{'log'}" )
          or die
          "ERROR : $sSubName :: Could not open $phCmdLineArgs->{'log'} file for writing.Reason : $!\n";
    }
    @aRequired = qw(input_file);
    foreach $sOption (@aRequired) {
        if ( !defined( $phCmdLineArgs->{$sOption} ) ) {
            printLogMsg( $ERROR,
                "ERROR : $sSubName :: Required option $sOption not passed" );
        }
    }
}

####################################################################################################################################################

# Description   : Used to handle logging of messages(errors and warnings) during the execution of the script
# Parameters    : level = can be ERROR, WARNING or INFO
#		  msg   = msg to be printed in the log file or to STDERR
# Returns       : NA
# Modifications :

sub printLogMsg {
    my ( $level, $msg ) = @_;
    if ( $level <= $DEBUG ) {
        print STDERR "$msg\n";
        print $logfh "$msg\n" if ( defined($logfh) );
        die "" if ( $level == $ERROR );
    }
}

__END__

#####################
# POD DOCUMENTATION #
#####################

=head1 NAME

fix_gene_symbols.pl - Script to remove upper case, duplicate gene symbols as well replace existing gene symbols in DB with new ones provided in the rules file

=head1 SYNOPSIS

# USAGE : perl fix_gene_symbols.pl -i <meta data file> [ -o <path to output dir>  -u <database username> -p <database password> -s <database server> -t <path to .tbl file> -l <path to log file> ]

	parameters in [] are optional

=head1 OPTIONS

	-i <input_file>	:	Path to input meta data file containing one row for a database. Mandatory
[	
	-o <output_dir>	:	Path to the output directory where the corrected tbl file will be created. Mandatory if change required in only tbl file

	-t <tbl_file>	:	Path to the original .tbl file to be corrected for gene symbols. Mandatory if change required in only tbl file

	-l <log>	: 	Path to log file. Optional

] 

=head1 DESCRIPTION



=head1 INPUT

	Required input is the meta data file with the following columns in a row for a database:
		[0] = Db name
		[1] = NCBI locus tag
		[2] = Path to rules file
		[3] = Path to gene symbols file
		[4] = Bioproject Id
		[5] = Organism name
		[6] = strain name
		[7] = Organism serotype
		[8] = Host
		[9] = Date
		[10]= Country
		[11]= Assembly method
		[12]= Genome Coverage
		[13]= Sequencing platform
		[14]= Contact person's name. Last name\sFirst name
		[15]= Contact person's email address
		[16]= Comma-separated author list. Each name should be Last name\sFirst name\sMiddle Initial
		[17]= Title of the publication
	For this script only column [0] and [3] needs to be have valid values.

	Provide the original tbl file dumped from manatee database. 

=head1 OUTPUT

	Corrected .tbl file 

=head1 AUTHOR

	Sonia Agrawal
	Senior Bioinformatics Software Engineer
	Institute for Genome Sciences
	School of Medicine
	University of Maryland
	Baltimore, MD - 21201
	sagrawal@som.umaryland.edu

==cut
