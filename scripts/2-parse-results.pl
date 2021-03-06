#!/usr/bin/perl

# Script for parsing results generated by 
# the AccessoryGenome.pl pipeline.
# By Jeff Gauthier - v1.0 - 7/2/2016

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

$align_dir = "./Results/Alignments";
system ("mkdir $align_dir");


### STEP 1 - CLEANUP AND CREATE OUTPUT FILES ######################

open (BATCH, "<", "batch.txt") or die $!; # batch file generated by main script
while (<BATCH>){
	chomp $_;
	@args = (split(" ", $_)); # store batch arguments in an array 
	system("mv ./Results/$args[0]_$args[1] $align_dir"); # clean-up the mess in Results folder
	unless (-e "./Results/$args[0].scores.csv") { # if output file has been created, no need to overwrite it again.
		system ("echo \"Protein\" >./Results/header.txt");				
		system ("cat ./Results/header.txt $align_dir/$args[0]_$args[1]/raw-list.txt > ./Results/$args[0].scores.csv");
		system ("paste ./Results/$args[0].scores.csv ./$args[0].seq/$args[0].functions >./Results/$args[0].temp.csv"); #append functions		
		system ("cat ./Results/$args[0].temp.csv >./Results/$args[0].scores.csv"); #overwrite file with temp
		system ("paste ./Results/$args[0].scores.csv ./$args[0].seq/$args[0].coords >./Results/$args[0].temp.csv"); #append coords		
		system ("cat ./Results/$args[0].temp.csv >./Results/$args[0].scores.csv"); #overwrite file with temp			
		system ("rm ./Results/header.txt");
		system ("rm ./Results/$args[0].temp.csv");	
	}
}
close BATCH;


### STEP 3 - APPEND IDENTITY SCORES TO OUTPUT FILE #################

open (BATCH, "<", "batch.txt") or die $!; # batch file generated by main script
while (<BATCH>) {
	chomp $_;
	@args = (split(" ", $_)); # store batch arguments in an array 
	system ("paste ./Results/$args[0].scores.csv $align_dir/$args[0]_$args[1]/results.txt > ./Results/$args[0].temp.csv");
	system ("cat ./Results/$args[0].temp.csv >./Results/$args[0].scores.csv");
	system ("rm ./Results/$args[0].temp.csv");
}
close BATCH;

### STEP 4 - FINISHING OUTPUT FILES AND QUIT #################

open (BATCH, "<", "batch.txt") or die $!; # batch file generated by main script
while (<BATCH>) {
	chomp $_;
	@args = (split(" ", $_)); # store batch arguments in an array
	system ("sed -i \"s/\.out//g\" ./Results/$args[0].scores.csv");
}
close BATCH;

print "Generating result files in ./Results/ .......Done!\n";
 
