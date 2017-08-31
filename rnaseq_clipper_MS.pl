#!/usr/bin/perl

$usage= "

rnaseq_clipper.pl  : 

Clips 5'-leader off Illumina fastq reads in RNA-seq
Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (mstudiva@fau.edu) 

Retains duplicated reads sharing the same degenerate header and removes 
the first 6 bases of the sequence (reads containing N bases in this
region are discarded, too)

prints to STDOUT

arguments:
1 : fastq file name
2 : string to define the leading sequence, default '[ATGC]{4}G{3}'
'keep' : optional flag to say whether the sequences without leader should be kept. 
		 By default, they are discarded.

Example:
rnaseq_clipper.pl D6.fq

					 
";

my $fq=shift or die $usage; # shift pulls the filename into variable fq, die kills script if no filename is given, then loops back to beginning of script
my $lead=""; # adaptor is undefined
my $keep=1; # keep is undefined, default is to remove sequences without headers
if ($ARGV[0]) { $lead=$ARGV[0];} # if a first argument is given, use it for the adaptor sequence
else { $lead="[ATGC]{4}G{3}";} # otherwise use the default
if ($ARGV[1]) {$keep=1;} # if second argument is given, sequences missing adaptor are kept
# specifying variable names to prevent autovivification
my $trim=0; # trim value set to undefined
my $name=""; # name is undefined
my $name2=""; # name 2 is undefined
my $seq=""; # sequence is undefined
my $qua=""; # quality score is undefined
my %seen={}; # seen is undefined
my $dups=0; # duplicates is undefined
my $tot=0; # total is undefined
my $nohead=0; # no headers is undefined
my $kept=0; # kept is undefined
open INP, $fq or die "cannot open file $fq\n"; # open input of given filename, or kills script with message
my $ll=3; # ll is initially set to three (option 1)
while (<INP>) { # while reading input files
	if ($ll==3 && $_=~/^(\@.+)$/ ) { # if option 1 and variable _ is something (.+) between the beginning (^) and end ($) of a string starting with @ (header of each read)
		$name2=$1;  # name2 is set to the read header starting with @
		$tot++; # returns the next consecutive value of total, aka counts total reads
		if ($seq=~/^($lead)(.+)/) {	# if the sequence starts with the adaptor, followed by something (read sequence)
			my $start=substr($2,0,6); # start variable is the first 6 characters of the read sequence
			my $idtag=$1.$start; # defines idtag as the adaptor.first 6 bases of sequence
			# if (!$seen{$idtag} and $idtag!~/N/) { # if the particular idtag is seen and does not include an N
			if ($idtag!~/N/) {
				if(!$seen{$idtag}){	
						$seen{$idtag}=1; # sets the count for that particular idtag equal to one
					}
				else{
						$seen{$idtag}=$seen{$idtag}+1;
				}
					$trim=length($1); # and trim is equal to the length of the adaptor
					$kept++; # adds +1 for every new adaptor seen
					print "$name\n$2\n+\n",substr($qua,$trim),"\n"; # prints read header, sequence, +, quality scores minus the corresponding adaptor regions, separated by new lines
					} 
			else { 
				$dups++; 
				} # or if non-unique adaptor seen, add +1 to duplicate count
		} 
		else { # if no adaptor is found
			$nohead++; # add +1 to nohead count
			if ($keep and $name) { print "$name\n$seq\n+\n",$qua,"\n";} # if keep and name are defined, prints name, sequence, +, quality scores with new lines in between
		}
		$seq=""; # sequence is undefined
		$ll=0; # ll is undefined
		$qua=""; # quality score is undefined
		@sites=(); # sites array is undefined
		$name=$name2; # sets the something string to be the original sequence
	}
	elsif ($ll==0){ # after the above, if ll is equal to zero (noheader, option 2)
		chomp; # remove any newline characters from the end of each read
		$seq=$_; # the sequence is kept the same without trimming (noheader kept)
		$ll=1; # ll set to one
	}
	elsif ($ll==2) { # or if ll is equal to two (option 3)
		chomp; # remove any newline characters from the end of each read
		$qua=$_; # quality scores kept the same without trimming (full sequence)
		$ll=3; # ll set to three
	}
	else { $ll=2;} # otherwise ll is set to two
}
# outputs the following strings with tabs in between:
# filename, total reads, missing header, duplicates, kept
warn " 
file\ttotal\tnoheader\tduplicates\tkept
$fq\t$tot\t$nohead\t$dups\t$kept

";
