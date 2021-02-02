#!/usr/bin/perl

#################################################################
# Accessory script to main script: AccessoryGenome.pl
# 1st step of the workflow: extracting protein sequences 
# from GenBank accessions and store in respective folders.
# By Jeff Gauthier - v1 - 7/2/2016
#--------------------------------------------------------
# Usage: ./scripts/1-gbk-extract.pl seq1.gbk seq2.gbk ... seqN.gbk
# Dependencies: EMBOSS v6.x.x.x (extractseq and transeq)
#
# INPUT VARIABLES
# Input file names are stored in @ARGV.
##################################################################

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


use warnings;

### CHECK IF EMBOSS IS HERE ############################################

system('embossversion 1> embv.txt 2>/dev/null');
open (EMBV, '<', 'embv.txt') or die $!;
if (<EMBV> eq "") {
	close EMBV;
	system ('rm embv.txt');
	die "EMBOSS not installed! Please install by typing \"sudo apt-get install emboss\" in your shell\n";
}
else {
	close EMBV;
	system ('rm embv.txt');
	print "EMBOSS installed on this machine.  Continuing...\n";
}


### EXTRACT CDS AND PROTEIN SEQUENCES #####################################
### Approach - pipeline for extractseq and transeq (EMBOSS)################

foreach my $gbk (@ARGV) {	#for each gbk file included as command-line arguments...

	print "Processing GenBank file $gbk...\n";	
	system("mkdir $gbk.seq");
	system ("mkdir $gbk.seq/CDS $gbk.seq/PEP");

	open (GENBANK, "<", "$gbk") or die $!;	#input GBK file
	open (FEATURES, ">>", "$gbk.seq/$gbk.features") or die $!;	#output features list
	
	while (<GENBANK>){	#read genbank file and keep only CDS coordinates.
		@lines = grep (/CDS  /, $_);
		print FEATURES join ("\n", @lines);	
	}
		
	close GENBANK;
	close FEATURES;

	# generate one-line genbank string for pattern matching
	open (GENBANK, "<", "$gbk") or die $!;
	open (TEMP, ">", "./$gbk.seq/$gbk.temp") or die $!;
	while (<GENBANK>){
		chomp $_;
		$temp .= $_;	
	}
	@temp_lines = $temp =~ /CDS.+?product=\".+?\"/g;
	print TEMP join ("\n", @temp_lines);
	close TEMP;
	close GENBANK;	

	#store CDS functions in a separate file
	open (TEMP, "<", "./$gbk.seq/$gbk.temp") or die $!;	
	open (FUNCTIONS, ">", "$gbk.seq/$gbk.functions") or die $!;
	print FUNCTIONS "Function\n";
	while (<TEMP>) {
		if (/CDS.+?product=\"(.+?)\"/) {print FUNCTIONS $1, "\n";}
	}	
	close TEMP;	
	close FUNCTIONS;

	#store coordinates of CDS in an array
	open (FEATURES, "<", "$gbk.seq/$gbk.features") or die $!;
	open (COORDS, ">", "$gbk.seq/$gbk.coords") or die $!; 
	while (<FEATURES>){	#remove tags and spaces and store in @features array 
		$_ =~ s/CDS//g;
		$_ =~ s/ //g;
		push (@features, $_);	
	}
	print COORDS "Coords\n"; #header	
	print COORDS join ("", @features);
	close COORDS;
	close FEATURES;

	$num = 1; #annotation numbering variable

	foreach my $i (@features) {
		
		if (grep (/complement/, $i)) {	#if the CDS needs reverse-complementation...
			chomp ($i);
			print "Extracting sequence $num ($i) of $gbk\n";
			#print "found $i no newline\n";
			$i =~ s/complement//;	#remove tags and brackets
			$i =~ s/join//;	#remove join comment
			$i =~ s/\(//g;
			$i =~ s/\)//g;
			#time to launch emboss commands
			system("extractseq -sequence ./$gbk -regions $i -outseq ./$gbk.seq/CDS/$gbk.$num.fna 2>/dev/null"); 
			system("revseq -sequence ./$gbk.seq/CDS/$gbk.$num.fna -outseq ./$gbk.seq/CDS/$gbk.$num.fna 2>/dev/null");
			system("transeq -sequence ./$gbk.seq/CDS/$gbk.$num.fna -table 11 -outseq $gbk.seq/PEP/$gbk.$num.tfa 2>/dev/null");
			$num = $num + 1;
		}
		else { # if it does not...
			chomp ($i);
			print "Extracting sequence $num ($i) of $gbk\n";
			#print "found $i no newline\n";
			$i =~ s/join//; #remove join comment
			$i =~ s/\(//g;
			$i =~ s/\)//g;
			#time to launch emboss commands
			system("extractseq -sequence ./$gbk -regions $i -outseq ./$gbk.seq/CDS/$gbk.$num.fna 2>/dev/null"); 
			system("transeq -sequence ./$gbk.seq/CDS/$gbk.$num.fna -table 11 -outseq $gbk.seq/PEP/$gbk.$num.tfa 2>/dev/null");
			$num = $num + 1;
		}
	}
	@features = (); #reset features array
}
$temp = "";
print "Done!\n";
