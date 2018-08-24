#!/usr/bin/perl
use warnings;
use strict;
use feature qw(say);
use Scalar::Util qw(looks_like_number);
use Getopt::Long;
use Pod::Usage;
use Data::Dumper qw(Dumper);
use File::Basename;



######################################################################################
#
# File	  : vcf2plink.pl
# History : 6/1/2018 Created by Kevin Freeman(KF)
#
#######################################################################################
#
# This script takes a vcf file and a corresponding population (.txt) file and combines
# and reformats the information, outputting a .ped and a .map file that can be used 
# in hapflk analysis or anything else requiring PLINK format
#
#######################################################################################

my ($vcf, $popFile); 
my ($outfile);

########### Command line options ######################################################

my $usage = "\nUsage: $0 [options]\n 
Options: 

     -vcf		VCF to convert  (required)
     
     -population	Corresponding population file (required)
     
     -outfile		Output file prefix for .ped and .map file (default: <name_of_vcf>)
     
     -help		Show this message

";

GetOptions(
   'vcf=s' => \$vcf,
   'population=s'	=> \$popFile,
   'outfile=s' 		=> \$outfile,
    help => sub { pod2usage($usage); },
) or pod2usage(2);

unless ($vcf) {
    die "\n-vcf not defined\n$usage";
}
unless ($popFile) {
	die "\n-pop not defined\n$usage";
}
unless (defined $outfile){
	$outfile =  basename($vcf, ".vcf");
}

############ Read pop file into a hash ################################################
say STDERR "\nReading pop file.....";

my $popFh;
unless (open ($popFh, "<", $popFile) ){
	die "can't open '$popFile', $!";
}

my %individualsHash;
my $individualNum = 1;				    # individualNum assigns an iid to each individual in order. Starts at 1 becaue
										# 0 is not a valid plink iid
my ($inFinal, $colPheno, $colGroup);
my @groups;
										
while(<$popFh>){
    chomp $_;
    my @line = split(" ",$_);
    unless (looks_like_number $line[0]) { # make sure the program doesn't store the header
    	($inFinal)  = grep { $line[$_] eq "\"infinal\"" }    0..$#line;
		($colGroup) = grep { $line[$_] eq "\"group\"" }      0..$#line;
		($colPheno) = grep { $line[$_] eq "\"phenotype.\"" } 0..$#line;
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
    push @groups, $group;							    # save all groups so they can be counted later
    
    if (int($group) != $group){
    	# make sure user gave the correct column number
    	die "Invalid group value: ",$group,", are you sure the headings match the data?";
    }
   
  	$individualsHash{$individualNum}{group}  = $group;
   
  	# add phenotypes if user gave column numbers, else phenotype = -9
  	if (defined $colPheno){
  		$individualsHash{$individualNum}{phenotype} = $line[$colPheno];
  	}
  	else {
  		$individualsHash{$individualNum}{phenotype} = -9;
  	}
   			
    $individualNum++;
}
# individualsHash is now a hash of hashes. The 'outer' hash has keys that correspond with
# the individual numbers and the inner, anonymous hash has keys that correspond with the type
# of data we want to access. It always has the key "group" to access all the groups of the individual.
# If the user defined them, it may also have keys for "phenotype" 

close $popFh;

############  Process VCF  ###################################################
say STDERR "Reading VCF.....";

my $vcfFh;
unless(open ($vcfFh, "<", $vcf)){
	die "Could not open $vcf";
}

my @snpValArray = vcf2ped($vcfFh);

close $vcfFh;
#### Convert vcf into map file
my $vcfFh2;
unless(open($vcfFh2, "<", $vcf)){
	die "Could not open $vcf $!";
}
my $mapOutFile = join(".", $outfile, "map");
vcf2map($vcfFh2,$mapOutFile);

close $vcfFh2;

# snpValArray is now an array of strings. Each string corresponds to one snp. It contains
# a list of alleles for each individual. Inviduals are separated by spaces and alleles
# within an individual are separated by commas 


########### Put data into .ped format ###############################################

say STDERR "Writing .ped file.....";

## create outfile name if one was not given, open outfile

my $pedOutFile = join(".", $outfile, "ped");

my $pedOutFh;
unless (open ($pedOutFh, ">", $pedOutFile)){
	die "Could not open $outfile $!";
}

# go through the hash in order, by individual
foreach my $individual (sort {$a <=> $b} keys %individualsHash ){
	my @allIndivAlleles;
	foreach my $snp (@snpValArray){  
		my $indivAllele  = $snp -> {$individual};	 	# access the value in the hash ref for the current individual
		$indivAllele =~ s/,/ /g;						# replace commas with spaces
		push @allIndivAlleles, $indivAllele;
	}

	say $pedOutFh join(" ", $individualsHash{$individual}{group}, $individual, 
					0, 0, 0, $individualsHash{$individual}{phenotype}, @allIndivAlleles);
}

close $pedOutFh;
say STDERR "Created $pedOutFile";
my @uniqueGroups = uniq(@groups);
say scalar @uniqueGroups;         		# ouput the number of groups in the dataset to stdout so it can be easily recorded programmatically   


################# SUBROUTINES ############################################################

#-----------------------------------------------------------------------
# void context = vcf2map($infh, $outfh);
#-----------------------------------------------------------------------
# This subroutine takes an infile handle and an outfile handle. It 
# takes the input vcf file, extracts the relevant data, and generates
# a .map file 
#-----------------------------------------------------------------------
sub vcf2map{
	my ($inHandle, $outFile) = @_;
	my $outFh;
	unless (open($outFh, ">", $outFile)){
		die "Could not open $outFile for writing";
	}
	my $i = 1;							# start at 1 because 0 is an invalid plink IID
	while (<$inHandle>){
		if ($_ =~ /^#/){				#skip the header lines
			next;
		}
		my @line  = split("\t", $_);
		my $chrom = $line[0];
		my $varID = $line[2];
		my $pos   = $line[1];
		
		if ($varID eq '.'){				# replace missing var id data with a unique number
			$varID = $i;
		}
		
		say $outFh join("\t", $chrom,$varID,"0",$pos); # print to the .map file
		$i++;
	}
	close $outFh;
	say STDERR "Created $outFile\n";
	
}

#-----------------------------------------------------------------------
# void context = vcf2ped($infh);
#-----------------------------------------------------------------------
# This subroutine takes an infile handle. It takes the input vcf file, 
# extracts the relevant data, and generates a line of alleles 
#-----------------------------------------------------------------------
sub vcf2ped {
	my ($infh) = @_;
	
	my @allelesConverted;
	
	while (<$infh>){
		if ($_ =~ /^#/){			#skip header lines
			next;
		}
		my @line = split("\t",$_);
			
		my $refAllele = $line[3];
		my $altAllele = $line[4];
	
		my $translatedLine = _vcfLine2basesLine( { reference => $refAllele, 
											  	   alt       => $altAllele, 
											  	   line 	     => $_ } );
		push @allelesConverted, $translatedLine;
	}
	return @allelesConverted;
}

#-----------------------------------------------------------------------
# ($commaDelimitedString) = _vcfLine2basesLine({reference, alt, line});
#-----------------------------------------------------------------------
# This subroutine takes a line containing genotype information and 
# parameters to convert that information in the form of a referenced
# hash and converts them to a line that represents the info as bases.
# It returns a reference to a hash where the keys are individuals 
# and values are the phenotypes of those individuals at the given SNP
#-----------------------------------------------------------------------

sub _vcfLine2basesLine {
	my ($refArgs) = @_;
	my $line = $refArgs 	-> {line};
	my $ref = $refArgs 		-> {reference};
	my $alt = $refArgs 		-> {alt};
	
	my @gtArray    = $line =~ /[\d\.]\|[\d\.]/g;     # get only the genotype data
	my $basesLine  = join(" ", @gtArray);		    		 	
 	
	$basesLine =~ s/0/$ref/eg;					     # replace 0's with ref allele
	$basesLine =~ s/1/$alt/eg;					     # replace 1's with alt allele
	$basesLine =~ s/\|/,/g;						     # replace |'s with commas
	$basesLine =~ s/\./?/g;						     # replace missing data with 0
	
	my @vals = split(" ", $basesLine);
	my @indivs = (1..scalar @vals);					 # create a sequence from 1 to the length of the array	

	my %snpHash;									
	@snpHash{@indivs} = @vals;						 # create a hash where the keys are the individuals and values are genotypes
	return \%snpHash;								 # return referenced hash
}
#-----------------------------------------------------------------------
# @uniqueArray = uniq(@arrayWithDups);
#-----------------------------------------------------------------------
# This subroutine removes any duplicate values from an array
#-----------------------------------------------------------------------
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}
