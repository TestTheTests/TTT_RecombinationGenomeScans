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
# File	  : splitPedMapbyChrom.pl
# History : 6/7/2018 Created by Kevin Freeman(KF)
#
#######################################################################################
#
# This script takes a .ped and a .map file and splits them into files that only contain
# information about one chromosome each 
#
#######################################################################################

my ($fileName); 
my ($outfile);

########### Command line options ######################################################

my $usage = "\nUsage: $0 [options]\n 
Options: 

	-infile		file prefix of .ped and .map files to be split 
	
	-outfile	Output file prefix for .ped and .map files (default: file name + chromosome number)
     
	-help		Show this message

";

GetOptions(
   'outfile=s' 	=> \$outfile,
   'infile=s'	=> \$fileName,
    help => sub { pod2usage($usage); },
) or pod2usage(2);

unless ($fileName) {
    die "\n You must enter the name of a set of files to convert";
}
unless (defined $outfile){
	$outfile = basename($fileName);
}

############ Read map file into a hash ################################################

say STDERR "\nReading .map file.....";

my $mapFh;
my $mapFile = join(".", $fileName, "map");
unless (open ($mapFh, "<", $mapFile) ){
	die "can't open $fileName, $!";
}

my %mapHash;
my $alleleNum = 0;

while(<$mapFh>){
    chomp $_;
    my @line = split("\t",$_);
    unless (looks_like_number $line[0]) { # make sure the program doesn't store the header
    	warn join("\t", "Format not recognized, Skipping line: @line");
    	next;
    }
 	my $chr  = $line[0];
 	my $id   = $line[1];
 	my $dist = $line[2];
 	my $pos  = $line[3];
 	
 	if (defined $mapHash{$chr}){								# if chr has been seen before, add info to arrays 
 		# get previous values for pos and index
 		my $indArrRef  = $mapHash{$chr}{index};
 		my $posArrRef  = $mapHash{$chr}{pos};
 		my $distArrRef = $mapHash{$chr}{dist};
 		my $idArrRef   = $mapHash{$chr}{id};
 		# push the new info
 		push @$indArrRef, ($alleleNum, $alleleNum + 1);
 		push @$posArrRef, $pos;
 		push @$idArrRef, $id;
 		push @$distArrRef, $dist; 
 		# add new arrays to hash
 		$mapHash{$chr}{index} = $indArrRef;
 		$mapHash{$chr}{pos}   = $posArrRef;
 		$mapHash{$chr}{id}	  = $idArrRef;
 		$mapHash{$chr}{dist}  = $distArrRef;
 	}
 	else {
 		$mapHash{$chr}{index} = [$alleleNum, $alleleNum + 1];
 		$mapHash{$chr}{pos}   = [$pos];
 		$mapHash{$chr}{dist}  = [$dist];
 		$mapHash{$chr}{id}	  = [$id];
 	}
 	$alleleNum += 2;
  
}
# %allelesHash is now a hash of hashes of array references. Chromosomes are used as the keys to access
# positions, indexes, and variant ids for each chromosome

close $mapFh;

# use the data structure with the map file information to print the new mapfiles
printMapFiles(\%mapHash, $outfile);

############ Read .ped file into an array ###################################################

say STDERR "\n\nReading .ped file.....";

my $pedFh;
my $pedFile = join(".", $fileName, "ped");
unless(open ($pedFh, "<", $pedFile)){
	die ("Could not open $pedFile");
}

chomp(my @pedArray = <$pedFh>);
close $pedFh;

# print new .ped files using @pedArray and %mapHash
printPedFiles(\@pedArray, \%mapHash, $outfile);

################# SUBROUTINES ############################################################

#-----------------------------------------------------------------------
# printMapFiles(\%mapHash, $outfile);
#-----------------------------------------------------------------------
# This subroutine takes a data structure representing the contents of a
# .map file in hashes and prints new map files, one for each chromosome 
#-----------------------------------------------------------------------
sub printMapFiles {
	my ($mapHashRef, $outfile) = @_;
	foreach my $chr (keys %$mapHashRef){
		# generate and open file of the form: <orig file name>_chr<>num>.map
		my $outChrFile = join("", $outfile, "_chr", $chr, ".map");
		my $outFh;
		unless (open($outFh, ">", $outChrFile)){
			die "Could not open $outChrFile";
		}
		# access info for current chromosome
		say "Creating $outChrFile";
		my $idsRef  = $mapHashRef -> {$chr}{id};
		my $posRef  = $mapHashRef -> {$chr}{pos};
		my $distRef = $mapHashRef -> {$chr}{dist};
		
		#print to file
		my $i = 0;
		foreach my $id (@$idsRef){
			say $outFh join("\t", $chr, $id, $distRef -> [$i], $posRef -> [$i]);
			$i++;
		}
		close $outFh;
	}
}

#-----------------------------------------------------------------------
# printPedFiles(\@pedArray, \%mapHash, $outfile);
#-----------------------------------------------------------------------
# This subroutine takes a data structure representing the contents of a
# .map file in hashes, an array representing the contents of a ped file,
# and an outfile prefix and creates new .ped files separated by chromosome 
#-----------------------------------------------------------------------
sub printPedFiles {
	my ($pedArrRef, $mapHashRef, $outfile) = @_;
	foreach my $chr (keys %$mapHashRef){
		# generate and open file of the form: <orig file name>_chr<>num>.ped
		my $outChrFile = join("", $outfile, "_chr", $chr, ".ped");
		my $outFh;
		unless (open($outFh, ">", $outChrFile)){
			die "Could not open $outChrFile";
		}
		# access info for current chromosome
		say "Creating $outChrFile";
		my $indexesRef = $mapHashRef -> {$chr}{index}; 
		my @indexesArr  = @$indexesRef;
		
		# loop through each individual
		foreach my $indLine(@$pedArrRef){
			my @splitLine   = split(/(\d) ([ATCG\?]) /, $indLine);
			my @bases = split(" ", $splitLine[3]);
			unshift @bases, $splitLine[2];                     # add base that we split on to array of other bases 
			my $lineStart = join("", @splitLine[0..1]);		   # add num that we split back onto the start of the line
			
			my @chrBases = @bases[@$indexesRef];			   # find all the bases that are on the current chr			
			
			# print to file
			say $outFh join(" ", $lineStart, @chrBases);
		}
		close $outFh; 
	}
}
