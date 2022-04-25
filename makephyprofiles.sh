#!/usr/bin/bash

# makephyprofiles.sh - convert OrthoFinder output to phylogenetic profiles
# Modified script from original by Dr. Daniel Barker, Computational and Evolutionary Genomics course 2022

# Process Orthofinder gene count output file to create tab-delimited
# phylogenetic profiles. This script generates one line of output per
# orthogroup. Each row contains the orthogroup name, followed by
# tab-delimited fields containing "1" if the orthogroup is present in
# the taxon and "0" otherwise.

# Run within the orthofinder directory
for iv in $(ls -d */); do
  trimmed=$(echo $iv | sed 's:/*$::')    # remove the trailing slash
	infile="Orthogroups.GeneCount.csv"     # input file (output by OrthoFinder)
  outfile="phyprofiles_${trimmed}.txt"   # specify output file
  echo -e "Orthogroup\tAcaryochloris_marina_strain_MBIC11017\tArthrospira_maxima_CS328\tCrocosphaera_watsonii_WH8501\tCyanothece_sp._Strain_PCC8801\tGloeobacter_violaceus_PCC7421\tLyngbya_sp._PCC8106\tMicrocystis_aeruginosa_strain_NIES843\tNodularia_spumigena_CCY9414\tNostoc_azollae_0708\tSynechocystis_sp._Strain_PCC6803" > $outfile             # insert header names
  # Ignore the first line, which contains column headings.
  # For subsequent lines, output phylogenetic profiles, which contain "0" for
  # a genome if the gene count is zero or "1" if the gene count is nonzero.
  # The final column in the Orthofinder gene count file is not relevant to
  # the task and is ignored.
	awk ' NR > 1 {output = $1; for (field = 2; field <= NF - 1; field++) {if ($field == 0) {output = output "\t" 0;} else {output = output "\t" 1;}} print output;} '\
   $trimmed/Results*/$infile >> $outfile
  echo -e "Processed $trimmed"
done

