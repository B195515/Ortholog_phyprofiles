#!/usr/bin/bash

# makephyprofiles.sh - convert OrthoFinder output to phylogenetic profiles

# Process Orthofinder gene count output file to create tab-delimited
# phylogenetic profiles. This script generates one line of output per
# orthogroup. Each row contains the orthogroup name, followed by
# tab-delimited fields containing "1" if the orthogroup is present in
# the taxon and "0" otherwise.

for iv in $(ls -d */); do
  trimmed=$(echo $iv | sed 's:/*$::');
	infile="Orthogroups.GeneCount.csv";      # input file (output by OrthoFinder)
  outfile="phyprofiles_${trimmed}.txt";
# Ignore the first line, which contains column headings.
# For subsequent lines, output phylogenetic profiles, which contain "0" for
# a genome if the gene count is zero or "1" if the gene count is nonzero.
# The final column in the Orthofinder gene count file is not relevant to
# the task and is ignored.
	awk ' NR > 1 {output = $1; for (field = 2; field <= NF - 1; field++) {if ($field == 0) {output = output "\t" 0;} else {output = output "\t" 1;}} print output;} '\
   $trimmed/Results*/$infile > $outfile;
done

