#!/usr/bin/perl -w
use strict;

## INPUT

# Get the input file (fasta or fastq) from first command line argument
chomp(my $inputFile = $ARGV[0]);

# Get the kmer size from second command line argument
chomp(my $kmerSize = $ARGV[1]);

# If the input file doesn not exist, ask if a default file should be used
my $defaultFile = "../data/cladonia_forward.fq";

unless (-e $inputFile) {
	print "Sorry, the file \"$inputFile\" does not exist.\n";
	# ask for using default file
	my $useDefaultFile = "";
	while($useDefaultFile ne "y" && $useDefaultFile ne "n") {
		print "Should I use the default file? [y/n]: ";
		chomp($useDefaultFile = <STDIN>);
	}
	if ($useDefaultFile eq "y") {
		$inputFile = $defaultFile;
	} elsif ($useDefaultFile eq "n") {
		print "Program exit.\n";
		exit;
	}
}


# subroutine to get sequences out of input file
sub getSequencesFromFile {

	# first argument is input file
	my $inputFile = $_[0];
	
	# Check if the input file is fasta or fastq, otherwise abort
	# therefor, look at the substring beginning after the last "." in the file name
	my $fileFormat = substr($inputFile, rindex($inputFile, '.')+1, length($inputFile));

	# array for saving sequences
	my @sequences;

	if ($fileFormat eq "fq" || $fileFormat eq "fastq" || $fileFormat eq "fa" || $fileFormat eq "fasta") {
		
		print "Collecting read sequences from \"$inputFile\"...\n";

		# get file contents
		open(IN, $inputFile)  or die "Unable to open file \"$inputFile\"";
		chomp(my @fileContents = <IN>);
		close(IN);

		# get sequences out of file contents and save thme into the array

		if ($fileFormat  eq "fq" || $fileFormat eq "fastq") {
			# in case of a fastq file, sequences are separated by 4 lines
			for(my $i = 1; $i < scalar(@fileContents); $i += 4) {
				push(@sequences, $fileContents[$i]);
			}
		} elsif ($fileFormat eq "fa" || $fileFormat eq "fasta") {
			# in case of a fasta file, sequences are separated by 2 lines
			for(my $i = 1; $i < scalar(@fileContents); $i += 2) {
				push(@sequences, $fileContents[$i]);
			}
		} else {
			print "The input file must be a fasta or fastq file\n";
			exit;
		}
	}
	
	# return the array with all sequences
	return @sequences;
}

# Hash table for kmer coverage
my %kmerHash;

# get all sequences from input file
my @sequences = getSequencesFromFile($inputFile);

print "Collecting k-mers from read sequences...\n";
# get all kmers ocurring in the sequences and fill in the kmer coverage hash table
foreach(@sequences) {
    
	# current read sequence in the iteration
	my $sequence = $_;
	
	# naive algorithm to get all kmers out of a sequence
	# iterates through the sequence by moving at each step one position to the right until last kmer reached
	
	# find k-mers of lenght $kmerSize in the read sequence
	for (my $i = 0; $i <= length($sequence)-$kmerSize; $i++) {
		my $currentKmer = substr($sequence, $i, $kmerSize);
		
		## IMPORTANT: ignore multiple ocurrences of a kmer within the same read sequence, counts just as one
		# check if the current kmer equals another kmer with a lower position  in the current read sequence
		# that means, the current kmer has already been registered
		# index returns the position of the first ocurrence in the read for the sequence of the current kmer
		# $i is the position of the current kmer in the read (current iteration index)
		if (index($sequence, $currentKmer) < $i) {
			# don't add the current kmer to the hash again
			# continue to next iteration step
			next;
		} else {
			# add kmer to hash table
			$kmerHash{"$currentKmer"} += 1;
		} 
	}
}


# Create output csv file containing the kmer ocurrences from hash table
my @hashKeys = keys(%kmerHash);
print "Ready to print k-mers hash table. Press any key to continue: ";
<STDIN>;
open(OUT, ">kmers.csv") or die "Unable to open file \"kmers.csv\"";
print OUT "'K-mer','Coverage'\n";
# print all kmers in hash table
for (my $i=0; $i < @hashKeys; $i++) { 
	print OUT "'".$hashKeys[$i]."',".$kmerHash{$hashKeys[$i]}."\n";
}
close(OUT);




