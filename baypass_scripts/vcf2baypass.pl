#!/usr/bin/perl
use warnings;
use strict;
use feature qw(say);
use Scalar::Util qw(looks_like_number);
use Getopt::Long;
use Pod::Usage;
use Data::Dumper qw(Dumper);

######################################################################################
#
# File	  : vcf2baypass.pl
# History : 2/2/2018 Created by Kevin Freeman(KF)
#			2/9/2018 KF added col_group functionality
#
#######################################################################################
#
# This script takes a vcf file and a corresponding population file, extracts the
# genotype information for each individual, groups the genotype informatiion,
#  and outputs a file in a format that can be used in BayPass
#
#######################################################################################

my ($vcf, $popFile); 
my ($inFinal, $colEnv, $colPheno, $colGroup);
my ($outGeno, $outCovar);

########### Command line options ######################################################

my $usage = "\nUsage: $0 [options]\n 
Options: 

     -vcf		VCF to convert  (required)
     
     -population	Corresponding population file (required)
     
     -colGroup		Specify the column of the population file that your 
     			group data is in. 1st col is 0 (required)
     			
     -colEnv		Specify the column in the population file that your
     			environment data is in (optional)
     			
     -colPheno		Specify the column in the population file that your
     			phenotype data is in (optional)
     -colInFinal	Specify the column in the population file that tells you whether
     			an individual is in the final dataset, ie your VCF (optional)
     			
     -outGeno		Output file prefix for genotype file (default: <name_of_vcf>)
     
     -outCovar		Output file prefix for covar file (default: <name_of_vcf>)
     
     -help		Show this message

";

GetOptions(
   'vcf=s' => \$vcf,
   'population=s'	=> \$popFile,
   'colGroup=i' 	=> \$colGroup,
   'colEnv=i'		=> \$colEnv,
   'colPheno=i'		=> \$colPheno,
   'colInFinal=i'	=> \$inFinal,
   'outGeno=s' 		=> \$outGeno,
   'outCovar=s'		=> \$outCovar,
    help => sub { pod2usage($usage); },
) or pod2usage(2);

unless ($vcf) {
    die "\n-vcf not defined\n$usage";
}
unless ($popFile) {
	die "\n-pop not defined\n$usage";
}
unless ($colGroup) {
	die "\n-col_group not defined\n$usage";
}

############ Read pop file into a hash ################################################
say STDERR "\nReading pop file.....";

my $popFh;
unless (open ($popFh, "<", $popFile) ){
	die "can't open '$popFile', $!";
}

my %groupsHash;
my $individualNum = 0;				    # individualNum tells us what line we are on. 0 = first line
										# AFTER the header
while(<$popFh>){
    chomp $_;
    my @line = split(" ",$_);
    unless (looks_like_number $line[0]) { # make sure the program doesn't store the header
    	next;
    }
    
    # if $inFinal was specified, check if the individual is in the final dataset. If not, skip
    if (defined $inFinal){
    	my $in = $line[$inFinal];
    	unless($in eq "TRUE"){
    		next;
    	}
    }
    
    my $group = $line[$colGroup];
    
    if (int($group) != $group){
    	# make sure user gave the correct column number
    	die "Invalid group value: ",$group,", are you sure you selected the right column?";
    }
    if (defined $groupsHash{$group}){		  # check if group matches a previously read group
   		$groupsHash{$group}{individuals}  = join(" ",$groupsHash{$group}{individuals}, $individualNum);
   		
   		# add environments and phenotypes if user gave column numbers
   		$groupsHash{$group}{environments} = join(" ",$groupsHash{$group}{environments}, $line[$colEnv])
   			if defined $colEnv;
   		$groupsHash{$group}{phenotypes}   = join(" ",$groupsHash{$group}{environments}, $line[$colEnv])
   			if defined $colPheno;
    }
    else {
    	$groupsHash{$group}{individuals} = $individualNum;
    	
    	# add environments and phenotypes if user gave column numbers
    	$groupsHash{$group}{environments} = $line[$colEnv]   if defined $colEnv;
    	$groupsHash{$group}{phenotypes}   = $line[$colPheno] if defined $colPheno;
    }
    $individualNum++;
}
# groupsHash is now a hash of hashes. The 'outer' hash has keys that correspond with
# the group numbers and the inner, anonymous hash has keys that correspond with the type
# of data we want to access. It always has the key "individuals" to access all the individual
# numbers. If the user defined them, it may also have keys for "environments" and "phenotypes"

close $popFh;

############ Read VCF into an array ###################################################
say STDERR "Reading VCF.....";

my $vcfFh;
unless (open ($vcfFh, "<", $vcf) ){
    die "cant open '$vcf', $!";
}

my @snpValArray;
while(<$vcfFh>){
	if ($_ =~ /[\d\.][\|\/][\d\.]/){				   # regex matches: 1|0 OR 1/0 OR .|. OR ./.
		my @snpVals = $_ =~ /[\d\.][\|\/][\d\.]/g;     # create array containing all matches
		my $snpValsString = join(" ", @snpVals); 
		push @snpValArray, $snpValsString;
	}
}
# snpValArray is now an array of strings. Each string is a space separated list of allele
# counts for one snp

close $vcfFh;

########### Put data into baypass format ###############################################

say STDERR "Converting.....";
my $baypass = "";
foreach my $snp (@snpValArray){
	my $alleles = calcAlleles($snp, \%groupsHash);
	$baypass = join("",$baypass,$alleles,"\n");
}

unless ($outGeno){
	$outGeno = $vcf;
}
unless ($outCovar){
	$outCovar = $vcf;
}
open (my $outFh, '>', $outGeno.".geno");
print $outFh $baypass;
close $outFh;

my $ngroups = keys(%groupsHash);
say STDERR "\nCreated ", $outGeno.".geno";

## create covariate file if env or pheno was defined
if (defined $colEnv or defined $colPheno){
	my $covarFile = $outCovar.".covar";
	printCovarFile(\%groupsHash, $covarFile);
	say STDERR "Created ", $covarFile, "\n\nVar 1 = environment, Var 2 = phenotype
(if both env and pheno are present)";
}

say STDERR "\nNumber of populations = ", $ngroups;
say $ngroups;

#-----------------------------------------------------------------------
# $alleleString = calcAlleles( $refArray, $refHash);
#-----------------------------------------------------------------------
# This function takes a referenced array and a referenced hash. The 
# array is ordered vcf data for one snp and the hash tells which 
# individual belongs to each group. Returns a string representing 
# population grouped counts for the snp
#-----------------------------------------------------------------------
sub calcAlleles{
	my ($snpVals, $groupsHashRef) = (@_);
	my %groupsHash = %$groupsHashRef;
	my @snpValsArray = split(" ", $snpVals);
	my $alleleString = "";
	unless (scalar @snpValsArray == $individualNum){
		die "VCF data does not match population data, different number of
		individuals", $!;
	}
	my @alleleArray;
	
	foreach my $group (sort { $a <=> $b } keys %groupsHash){
		my @individuals = split (" ", $groupsHash{$group}{individuals}); # find all individuals in the group
		my $groupVals = join(" ",@snpValsArray[@individuals]);	         # find values by index
		my ($allele1, $allele2) = countAlleles($groupVals);
		push @alleleArray, ($allele1,$allele2);
	}
	return join(" ", @alleleArray);
}

#-----------------------------------------------------------------------
# ($count1, $count2) = countAlleles(@values);
#-----------------------------------------------------------------------
# This function takes a string with space delimited vcf data in the form 
# "0|1" and it counts the total number of alleles  of each type
#-----------------------------------------------------------------------
sub countAlleles{
	my ($groupVals) = @_;
	my $count1 = () = $groupVals =~ /0/g;		# find number of matches of '0'
	my $count2 = () = $groupVals =~ /1/g;
	return ($count1, $count2);
}

#-----------------------------------------------------------------------
# printCovarFile(\%groupHash, $fileName;
#-----------------------------------------------------------------------
# This function, called in void context, takes a hash reference
# with covariate data and a file name and prints the covariate data to 
# the named file.
#-----------------------------------------------------------------------
sub printCovarFile{
	my ($hashRef, $fileName) = @_;
	my $covarData ="";
	my %covarHash = %$hashRef;
	
	my @sortedKeys = sort {$a <=> $b} keys %covarHash;
	
	if (defined $colEnv){
		foreach my $group (@sortedKeys){ # go through each group and get the average Env		
			my $environments      = $covarHash{$group}{environments};
			my @environmentsArray = split(" ", $environments);
			$covarData = join(" ", $covarData, getAverage(@environmentsArray));
		}
		$covarData = join("",$covarData, "\n"); # add a new line to separate covariates
	}
	if (defined $colPheno){
		foreach my $group (@sortedKeys){    # go through each group and get the average pheno values
			my $phenotypes		  = $covarHash{$group}{phenotypes};
			my @phenotypesArray   = split(" ", $phenotypes);
			$covarData = join(" ", $covarData, getAverage(@phenotypesArray));
		}
	}
	
	$covarData =~ s/^\s+//mg;					# trim leading whitespaces
	# print the data to the file
	open (my $covFh, '>', $fileName);
	say $covFh $covarData;
	close $covFh;
}

#-----------------------------------------------------------------------
# $average = getAverage(@values);
#-----------------------------------------------------------------------
# This function takes an array of numbers and returns the average in 
# scalar context 
#-----------------------------------------------------------------------
sub getAverage {
	my @numbers = @_;
	my $length = @numbers;
	my $sum = 0;
	
	foreach my $number(@numbers){
		$sum += $number;
	}
	return $sum/$length;
}
