#!/usr/bin/perl
use warnings;
use strict;
use feature qw(say);
use List::MoreUtils 'first_index';
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

######################################################################################
#
# File	  : pruneSnps.pl
# History : 4/6/2018  Created by Kevin Freeman(KF)
#	  : 4/11/2018 KF changed logic to parse header instead of using user input
#	  : 11/5/2018 no longer relies on position column, but now requires both vcf
#	              and scan_results to be in the same order
#######################################################################################
#
# This script takes a vcf file and a file containing scan results and returns a vcf
# with only the quasi-independent alleles. 
#
#######################################################################################

my $vcf;
my $scanResults;
my $outFile;

########### Command line options ######################################################

my $usage = "\nUsage: $0 [options]\n 
Options:
     -vcf		VCF to prune  (required)
     -scan		File with scan results (required)
     -outfile		Name of output file (default: <name_of_vcf>.pruned.vcf)
     -help		Show this message

";

GetOptions(
   'vcf=s' 			=> \$vcf,
   'scan=s' 		=> \$scanResults,
   'outfile=s' 		=> \$outFile,
    help => sub { pod2usage($usage); },
) or pod2usage(2);

unless ($vcf) {
    die "\n-vcf not defined\n$usage";
}
unless ($scanResults) {
	die "\n-scan not defined\n$usage";
}

######## Read in Scan Results ####################################################
my $scanFh;
unless (open($scanFh, '<', $scanResults)){
	die "Could not open $scanResults for reading $!";
}
my $first = 1;
my @indArray  = [];
my $indCol;
my $keepCol;
while (<$scanFh>){
	chomp $_;
	my @line = split(/\s/, $_); #split on whitespace
	if ($first){    # process header line
		$indCol    = first_index { /quasi_indep/ } @line; 
		$keepCol   = first_index { /keep_loci/ } @line;
		$first = 0;
	}
	else {
		my $ind  = $line[$indCol];
		my $keep = $line[$keepCol];
		
		if ($keep eq "TRUE"){
			push @indArray, $ind;
		}
	}
}

###### Read in VCF, check lines, print ###########################################

my $noExt = $vcf;
$noExt =~ s/\.vcf(\.gz|)//;    

say $noExt;
unless (defined $outFile){
	$outFile = join(".", $noExt, "pruned", "vcf");
}

# check if vcf is gzipped, if it is unzip it
my $inVCFfh;
if ($vcf =~/\.gz$/ ){
	`gunzip $vcf`;
	$vcf = join(".", $noExt, "vcf");
}
# open the files
unless (open($inVCFfh, '<', $vcf)){
	die "Can't open $vcf file for reading $!";
}
my $outFh;
unless (open($outFh, '>', $outFile)){
	die "Can't open $outFile for writing $!";
}

my $i = 0;
while (<$inVCFfh>){
	chomp $_;
	if ($_ =~ /^#/){        # print all header lines automatically	
		 say $outFh $_;   
		 next;
	} 
	my @line = split("\t", $_);
	
	if ($indArray[$i] eq "TRUE"){
		say $outFh $_;
	}
	$i++;
}

close $outFh;
close $inVCFfh;
`gzip $vcf`;
`gzip $outFile`;
my $nlines = $i + 1;

# if the vcf and scan results are not the same the pruned file is totally invalid. check lengths and delete the file if necessary
unless($nlines == scalar @indArray){
	`rm $outFile.gz`;
	die join(" ", "Length of scan results (", scalar @indArray, ") does not match length of vcf (", $nlines , "), $!"); 
}

say "\nCreated $outFile.gz";
