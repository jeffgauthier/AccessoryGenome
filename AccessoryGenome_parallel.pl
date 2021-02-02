#!/usr/bin/perl

###########################################################
# Script that finds shared and unique CDS genes among pairs of genomes.
# By Jeff Gauthier - v1.0 - 7/2/2016
# Usage: ./AccessoryGenome_parallel.pl -T <threads> genome1.gbk genome2.gbk ... genomeN.gbk
# Dependencies: EMBOSS and tfastx from FASTA v36.X.X suite.
###########################################################

###########################################################
# VERSION HISTORY
# v1.0-parallel - Homology search parallelized using Parallel::Loops module.
# v1.0 - Fully operational pipeline w/ annotations in result files.
# v0.9999 - Simplified pipeline.  No annotation feature for now.
###########################################################


###########################################################
# LICENSE

# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#############################################################

# COMMAND-LINE ARGUMENTS
# Stored in @ARGV.

use warnings;
use Math::Combinatorics;
use Getopt::Long;

GetOptions(
	"T=i" => \(my $maxProcs = 1),
) or die "________________________________________\nAccessoryGenome_parallel.pl v1.0\nScript that finds shared and unique genes among pairs of genomes.\nUsage: ./AccessoryGenome_parallel.pl -T <threads> genome1.gbk genome2.gbk ... genomeN.gbk\nDependencies: EMBOSS and tfastx from FASTA v36.X.X suite.\n\n";


# FASTA TOOL USED (tfastx36 ~ faster) (tfasty36 ~ slower but corrects frameshifts -- useful for 454 WGS contigs)
$engine = "tfastx36";

# E-VALUE THRESHOLD FOR REMOVING FALSE POSITIVE ALIGNMENTS
$threshold = "1e-10";

### PARALLELIZATION OF LOOPS ###############################

use Parallel::ForkManager;
$pm = new Parallel::ForkManager($maxProcs);


### MAIN CODE BLOCK ########################################

# STEP ZERO - Are there more than 2 genomes to compare?
unless ($#ARGV >= 1) {die "________________________________________\nAccessoryGenome_parallel.pl v1.0\nScript that finds shared and unique genes among pairs of genomes.\nUsage: ./AccessoryGenome_parallel.pl -T <threads> genome1.gbk genome2.gbk ... genomeN.gbk\nDependencies: EMBOSS and tfastx from FASTA v36.X.X suite.\n\n"
}

#startup message
print "________________________________________\nAccessoryGenome_parallel.pl v1.0\nScript that finds shared and unique genes among pairs of genomes.\nUsage: ./AccessoryGenome_parallel.pl -T <threads> genome1.gbk genome2.gbk ... genomeN.gbk\nDependencies: EMBOSS and tfastx from FASTA v36.X.X suite.\n\n";


#1st step - extracting CDS and Protein sequences
foreach my $gbk (@ARGV) {
	system ("./scripts/1-gbk-extract.pl $gbk");
	system("seqret -sequence $gbk -osf fasta -outseq ./$gbk.seq/$gbk.fasta 2>/dev/null"); # generate library file
	system("ls -v ./$gbk.seq/PEP > ./$gbk.seq/$gbk.pep.list"); #generate queries list
	}

# 2nd step - making a list of all possible pairwise genome searches
my $set = Math::Combinatorics -> new (count => 2, data => [@ARGV]); # combinatorics module
open (BATCH, ">", "batch.txt") or die $!;
while (my @combo = $set->next_combination) { # print permutations of @ARGV elements in BATCH file
	print BATCH "$combo[0] $combo[1]"."\n";
	print BATCH "$combo[1] $combo[0]"."\n";
	}
close BATCH;


# 3rd step - store BATCH elements in arrays and execute fasty36
system("mkdir ./Results");
open (BATCH, "<", "batch.txt") or die $!;
while (<BATCH>) {

	#store fasty inputs from each BATCH line in array
	chomp $_;
	@args = (split(" ", $_));

	print "Searching homologies in ", join (" x ", @args), "\n";	
	system("mkdir ./Results/$args[0]_$args[1]");
	system("mkdir ./Results/$args[0]_$args[1]/raw");
	#execute tfasty36 using arguments stored in array
	open (QUERIES,"<","./$args[0].seq/$args[0].pep.list");
	while (<QUERIES>) {
		my $pid = $pm->start and next;
		chomp $_;
		system("$engine -b 1 -d 1 -E $threshold -s BP62 -T 1 ./$args[0].seq/PEP/$_ ./$args[1].seq/$args[1].fasta > ./Results/$args[0]_$args[1]/raw/$_.out 2>/dev/null");		
		print "...";
		$pm->finish;		
	};
	$pm->wait_all_children;	
	close QUERIES;
	print "done!\n";

	# 4th step - parse TFASTY results

	#make list of results to parse in current folder	
	$folder = "./Results/$args[0]_$args[1]";
	system ("ls -v ./$folder/raw >> ./$folder/raw-list.txt");

	# open each file in list and store Results in a hash
	# keys = CDS name, values = Identity percentage
	open (LIST, "<", "./$folder/raw-list.txt") or die $!;
	open (RESULTS, ">", "./$folder/results.txt");
	print RESULTS "$args[1]_\%id\n"; #header
	while (<LIST>) {
		chomp $_; #remove newline
		open (FILE, "<", "./$folder/raw/$_");
		while (<FILE>) {
			if (/(\d+.\d)% identity/o){ 				
				print RESULTS "$1", "\n"; #print matching pattern			
				last;
			}
			elsif (/No sequences with E/){
				print RESULTS "N/S", "\n"; #print (N)ot/(S)ignificant
				last;
			}		
		}
	}
		close FILE;
	close LIST;
	close RESULTS;
	}

close BATCH;

### FINAL STEP -- SUMMARIZING RESULTS ####################

system("./scripts/2-parse-results.pl");
print "Run is over!  Back to shell...\n";
system("rm batch.txt");
