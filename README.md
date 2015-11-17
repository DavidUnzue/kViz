# kViz (k-mer coverage visualization)
kViz is a simple bioinformatics tool for visualizing the dristibution of the k-mer coverage within a given set of sequencing reads.

## Pre-Requisites
Since kViz is written in Perl, you'll need Perl 5.6.0 or later installed in your machine.

Additionally you'll need the [Template Toolkit](http://www.template-toolkit.com/) module for Perl. You can install it using the CPAN Module by just typing the following in the command line:
```
$ sudo cpan Template
```

## Usage
Run in the command line the Perl script under the *src* folder with all three arguments.
```
$ perl src/kviz.pl <INPUT_FILE> <OUTPUT_FILE_PREFIX> <K-MER_SIZE>
```

### Arguments
<INPUT_FILE>
    input fasta or fastq file containing read sequences
<OUTPUT_FILE_PREFIX>
    the path to the output file without file format ending, e.g.: "data/kmer_coverage_file"
<K-MER_SIZE>
    size of the created k-mers

## Output
kViz generates a HTML file with a histogram showing the k-mer coverge distribution. The histogram is generated using [Google Charts Tools](https://developers.google.com/chart/).

## Copyright
kViz is licensed under the MIT License, see [LICENSE.md](LICENSE.md) for more details.
