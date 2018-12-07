#!/usr/bin/perl
use warnings;
use strict;
use feature qw(say);
use Scalar::Util qw(looks_like_number);
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

######################################################################################
#
# File	  : getTableOfBaypassResults.pl
# History : 4/18/2018  Created by Kevin Freeman(KF)
#
#######################################################################################
#
# This script is used in the bash script 'run_baypass.sh' to compile the results of 
# all the analyses into a table
#
#######################################################################################

my ($start, $finish, $indir, $outdir);

########### Command line options ######################################################

my $usage = "\nUsage: $0 [options]\n 
Options:
     -start		Number of first VCF in batch (ex 10900)
     -finish		Number of final VCF in batch (ex 10960)
     -in                Where the vcf is stored     (default: ../../../results_final)
     -out               Where to write the results  (default: ../../../results)
     -help		Show this message

";

GetOptions(
   'start=i' 		=> \$start,
   'finish=i' 		=> \$finish,
   'in=s'		=> \$indir,
   'out=s'		=> \$outdir,
    help => sub { pod2usage($usage); },
) or pod2usage(2);

unless ($start) {
    die "\n-start not defined\n$usage";
}
unless ($finish) {
	die "\n-end not defined\n$usage";
}
unless ($outdir) {
	$outdir = "../../../results/";
}
unless ($indir) {
	$indir  = "../../../results_final/"; 
}
`mkdir -p $outdir`;

########## Loop through all results #######################################################
for (my $i = $start; $i <= $finish; $start++){
	say "\n$i\n";
	my $vcfFile = $indir."/".$i."_Invers_VCFallFILT.vcf";
	`gunzip $vcfFile.gz`; 
	unless (-e $vcfFile){ # make sure file exists
		say $i."_Invers_VCFallFILT.vcf not found, skipping";
		$i++;
		next;
	}
	open(my $outFh, '>', $outdir."/".$i."_baypass_results.txt")
		or die "Could not open outfile for $i";
	my $columnsHashRef = getColumns($i);
	printResults($columnsHashRef, $outFh);
	close $outFh;
	say "Created ".$outdir."/".$i."_baypass_results.txt";	
	`gzip $vcfFile`;    # recompress vcf for space 
	$i++;
}

#################### SUBROUTINES #############################################################

#-----------------------------------------------------------------------
# %colsHash = getColumns($number);
#-----------------------------------------------------------------------
# This function takes the number corresponding to a dataset and finds 
# the statistics of interest, which are spread across a number of files. 
# It returns these statistics as an hash of array refs, where the 
# referenced arrays represent columns of data 
#-----------------------------------------------------------------------

sub getColumns{
	my ($number) = @_;
	## find position
	my @posCol;
	my $vcfFile = $indir."/".$number."_Invers_VCFallFILT.vcf";
	if (open(my $vcfFh, '<', $vcfFile)){
		while (<$vcfFh>){
			my $line = $_;
			if ($line =~ /^##|POS/){
				next;
			}
			else {
				my @splitLine = split("\t", $line);
				push @posCol, $splitLine[1];
			}
		}
	}
	else {
		warn "Could not read pos from $vcfFile $!";
	}
	## find BF cols (stored as array refs)
	my ($ALLEnvColRef, $ALLPhenoColRef)  	   = _readBFfile('', $number);
	my ($PRUNEDEnvColRef, $PRUNEDPhenoColRef)  = _readBFfile("_ALL_PRUNED_MAT", $number);
	## find XtX cols (stored as array refs)
	my $ALLxtxColRef    = _readXtXfile('', $number);
	my $PRUNEDxtxColRef = _readXtXfile('_ALL_PRUNED_MAT', $number);
	
	# return all columns as references in a hash
	return { pos       => \@posCol, 
		 	 ALLEnv    => $ALLEnvColRef,    ALLPheno    => $ALLPhenoColRef,
		 	 PRUNEDEnv => $PRUNEDEnvColRef, PRUNEDPheno => $PRUNEDPhenoColRef,
		 	 ALLxtx    => $ALLxtxColRef,    PRUNEDXtx   => $PRUNEDxtxColRef, 
	};

}

#-----------------------------------------------------------------------
# my ($pheno, $env) = _readBFfile($prefix, $number);
#-----------------------------------------------------------------------
# This function takes a prefix and a number and parses the corresponding 
# BF file  to find the environment and phenotype BF values. It returns 
# two array references 
#-----------------------------------------------------------------------
sub _readBFfile{
	my ($prefix, $number) = @_;
	my @EnvCol;
	my @PhenoCol;
	my $bfFile = $number.$prefix."_summary_betai_reg.out";
	if (open (my $bfFh, '<', $bfFile)){
		while (<$bfFh>){
			my $line = $_;
			$line =~ s/^\s+//;							# trim leading whitespace
			my @splitLine = split(/\s+/, $line);		# split on 1 or more whitespace chars
			unless (looks_like_number($splitLine[0])){
				next; # skip header
			}
			my $BF = $splitLine[3];
			if ($splitLine[0] eq '1'){
				push @EnvCol, $BF;
			}
			elsif ($splitLine[0] eq '2'){
				push @PhenoCol, $BF;
			}
			else {
				warn "Unknown covariable number: $splitLine[0]";
			}
		}
	}
	else {
		warn "Could not read BF from $bfFile";
	}
	return (\@EnvCol, \@PhenoCol);
}

#-----------------------------------------------------------------------
# my ($xtxCol) = _readXtXfile($prefix, $number);
#-----------------------------------------------------------------------
# This function takes a prefix and a number and parses the corresponding 
# XtX file  to find the xtx values. It returns one array ref  
#-----------------------------------------------------------------------
sub _readXtXfile {
	my ($prefix, $number) = @_;
	my @xtxCol;
	if (open(my $xtxFh, '<', $number.$prefix."_summary_pi_xtx.out")){
		while (<$xtxFh>) {
			my $line = $_;
			$line =~ s/^\s+//;					    # trim leading whitespace
			my @splitLine = split(/\s+/, $line);	# split on whitespace
			unless (looks_like_number($splitLine[0])){ # skip header
				next;
			}
			push @xtxCol, $splitLine[5];		
		}
	}
	else {
		warn "Could not read xtx from summary for $number $!";
	}
	return \@xtxCol;
}

#-----------------------------------------------------------------------
# printResults(%resultsHash, $outFh);
#-----------------------------------------------------------------------
# This function takes the all the results from the different files in 
# the form of a hash and an open file handle to print to and prints the
# information to the file. It is called in void context
#-----------------------------------------------------------------------
sub printResults{
	my ($resultsHashRef, $outFh) = @_;
	say $outFh join("\t", "position","baypass_2.1_ALL_BF_pheno","baypass_2.1_ALL_BF_env",
					"baypass_2.1_ALL_XTX", "baypass_2.1_PRUNED_BF_pheno",
					"baypass_2.1_PRUNED_BF_env", "baypass_2.1_PRUNED_XTX");
	my $posRef         = $resultsHashRef -> {pos};
	my $allPhenoRef    = $resultsHashRef -> {ALLPheno};
	my $allEnvRef      = $resultsHashRef -> {ALLEnv};
	my $allxtxRef      = $resultsHashRef -> {ALLxtx};
	my $prunedPhenoRef = $resultsHashRef -> {PRUNEDPheno};
	my $prunedEnvRef   = $resultsHashRef -> {PRUNEDEnv};
	my $prunedXtxRef   = $resultsHashRef -> {PRUNEDXtx};
	
	my $i = 0;
	foreach my $pos( @$posRef){
		say $outFh join ("\t", $pos, @$allPhenoRef[$i], 
				  		 @$allEnvRef[$i], @$allxtxRef[$i], 
				  		 @$prunedPhenoRef[$i], @$prunedEnvRef[$i], 
				  		 @$prunedXtxRef[$i]);
		$i++;
	}
	
}
