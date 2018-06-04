#!/usr/bin/perl
use warnings;
use strict;
use feature qw(say);
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

######################################################################################
#
# File	  : pruneSnps.pl
# History : 4/6/2018  Created by Kevin Freeman(KF)
#		  : 4/11/2018 KF changed logic to parse header instead of using user input
#
#######################################################################################
#
# This script takes a vcf file and a file containing scan results and returns a vcf
# with only the quasi-independent alleles. It also generates a tab-delimited file 
# giving the indexes of the quasi-independent alleles
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
my %posHash;
my ($posCol, $chrCol, $indCol);
while (<$scanFh>){
	chomp $_;
	my @line = split(/\s/, $_); #split on whitespace
	if ($first){    # process header line
		($posCol) = grep { $line[$_] eq "\"pos\"" } 0..$#line;
		($chrCol) = grep { $line[$_] eq "\"chrom\"" } 0..$#line;
		($indCol) = grep { $line[$_] eq "\"quasi_indep\"" } 0..$#line;
		$first = 0;
	}
	else {
		my ($chr,$pos,$ind) = @line[$chrCol, $posCol, $indCol];
		# combine chr and pos into one key for the hash that maps to the 
		# independence 
		# ex : 1.37685 => "TRUE"
		#     chr pos       ind
		$chr =~ s/"//g;		# remove quotation marks around chr
		my $key = join(".", $chr,$pos); 
		$posHash{$key} = $ind;
	}
}

###### Read in VCF, check lines, print ###########################################

unless (defined $outFile){
	my $noExt = $vcf;
	$noExt =~ s/\.vcf$//;    # remove .vcf from the end of the vcf file name
	$outFile = join(".", $noExt, "pruned", "vcf");
}

my $inVCFfh;
unless (open($inVCFfh, '<', $vcf)){
	die "Can't open $vcf file for reading $!";
}
my $outFh;
unless (open($outFh, '>', $outFile)){
	die "Can't open $outFile for writing $!";
}
my $indexesOutFh;
my $indexFile = "indexes_remaining.txt";
unless (open ($indexesOutFh, '>', $indexFile)){
	die "Can't open indexes outfile for writing $!";
}

my $i = 0;
while (<$inVCFfh>){
	chomp $_;
	if ($_ =~ /^#/){        # print all header lines automatically	
		 say $outFh $_;   
		 next;
	} 
	my @line = split("\t", $_);
	
	my ($chr, $pos) = @line[0,1];
	my $key = join(".", $chr, $pos);
	
	if (defined $posHash{$key} and $posHash{$key} eq "TRUE"){
		say $outFh $_;
		print $indexesOutFh $i."\t"; 
	}
	$i++;
}

close $outFh;
close $inVCFfh;

say "\nCreated $outFile";
say "Created $indexFile\n";
