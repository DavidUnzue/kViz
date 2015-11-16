#!/usr/bin/perl -w
use strict;
# template toolkit
use Template;
# module for getting path info of a file
use File::Basename;

## INPUT

# print help if no arguments or just "-h" given
if (scalar(@ARGV) == 0 || $ARGV[0] eq "-h") {
	print "
	USAGE:
		perl kviz <INPUT_FILE> <OUTPUT_FILE_PREFIX> <K-MER_SIZE>

	ARGUMENTS:
		<INPUT_FILE> => input fasta or fastq file containing read sequences
		<OUTPUT_FILE_PREFIX> => the path to the output file without file format ending, format will be .csv
		<K-MER_SIZE> => size of the created k-mers
	","\n";
	exit;
}

# Get the input file (fasta or fastq) from first command line argument
chomp(my $inputFile = $ARGV[0]);

# Get the output file prefix from second command line argument
chomp(my $outputFilePrefix = $ARGV[1]);

# Get the kmer size from third command line argument
chomp(my $kmerSize = $ARGV[2]);

## BODY

# exit in case input does not exist
unless (-e $inputFile) {
	print "Sorry, the file \"$inputFile\" does not exist.\n";
	print "Program exit.\n";
	exit;
}

# create kmer abundance hash for input file
my %kmerHash = createKmerCoverageHash(getSequencesFromFile($inputFile));

# prepare data to be sent to template
my %templateData = ( col1 => 'kmer',
			col2 => 'coverage',
			kmers => \%kmerHash );
my $templateFilePath = dirname(__FILE__);
# initialize template object
my $tt = Template->new({
	INCLUDE_PATH => $templateFilePath.'/templates', #look here for template files
	RELATIVE => 1 #INCLUDE_PATH is relative
});
# send data to template object and create output file
$tt->process('histogram_template.tt2', \%templateData, $outputFilePrefix.'.html') || die $tt->error;


## SUBROUTINES

# subroutine to get sequences out of input file
sub getSequencesFromFile {

	# first argument is input file
	my ($inputFile) = @_;

	# Check if the input file is fasta or fastq, otherwise abort
	# therefor, look at the substring beginning after the last "." in the file name
	my $fileFormat = substr($inputFile, rindex($inputFile, '.')+1, length($inputFile));

	# array for saving sequences
	my @sequences = ();

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


# sub routine to create kmer coverage from a list of read sequences
sub createKmerCoverageHash {

	# first argument is an array with sequences
	my (@sequences) = @_;

	# hash table to write kamers and corresponding coverage to
	my %kmerHash = ();

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

	# return the hash with kmers and coverage
	return %kmerHash;
}
