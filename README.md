ACCESSORYGENOME.PL

Script that finds shared and unique CDS genes among pairs of genomes.
By Jeff Gauthier - v1.0 - 7/2/2016

Usage: ./AccessoryGenome_parallel.pl -T <threads> genome1.gbk genome2.gbk ... genomeN.gbk

Dependencies: 
 - EMBOSS
 - tfastx36 from FASTA v36.X.X suite
 - Perl modules Parallel::ForkManager and Math::Combinatorics


VERSION HISTORY

v1.0-parallel - Parallelized homolog searches (use -T parameter)
v1.0 - Fully operational pipeline w/ annotations in result files.
v0.9999 - Simplified pipeline.  No annotation feature for now.


LICENSE

This program is free software: you can redistribute it and/or modify  
it under the terms of the GNU General Public License as published by  
the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
General Public License for more details.

You should have received a copy of the GNU General Public License 
along with this program. If not, see <http://www.gnu.org/licenses/>.

