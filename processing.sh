# count number of unique sequences per file
for name in $(ls); do
  echo $name; 
  grep ">" $name | wc -l; 
done

grep -ni "^>" *.fasta | grep -ni "nif" > grepNif.txt

# get unique ids of proteins per species
for name in $(ls *.fasta); do
  echo $name >> fastasummary.txt;
  grep ">" $name | cut -c 2-9 | sort | uniq >> fastasummary.txt;
  echo -e "\n" >> fastasummary.txt;
done

# Run orthofinder for various inflation values
for value in {1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,4,6,10,12}; do
  # Run while in orthofinder folder
  dir="iv${value}"
  if [ -d "$dir" ]; then
    echo "$dir FOUND"
    continue
  fi
  mkdir iv$value
  cp ../proteinsets/*.fasta iv$value
  orthofinder -f iv$value -t 20 -I $value -og
done

# compare Statistics_Overall
echo "" > NumberofGenes.txt
for iv in $(ls); do
  echo $iv >> NumberofGenes.txt;
  head -n17 $iv/Results*/Statistics_Overall.csv >> NumberofGenes.txt;
  echo -e "\n" >> NumberofGenes.txt;
done

# search in Orthogroups.csv for nifH / dinitrogenase
echo "" > nifsearch.txt
for iv in $(ls); do
  echo $iv >> nifsearch.txt;
  grep -ni 'nif' $iv/Results*/Orthogroups.csv >> nifsearch.txt;
  echo -e "\n" >> nifsearch.txt;
done



# combine all 10 protein sets
for file in $(ls *.fasta); do
  cat $file >> allproteins.fasta
done
# manually remove any duplicate IDs
grep ">" allproteins.fasta | sort | uniq -c > file.txt # count repetitive IDs
grep -v "1 >" file.txt # show IDs with count > 1
grep -n "<ID>" allproteins.fasta # show which line has the repetitive IDs
sed 'm,n!d' file # print the lines from m to n
sed -i 'm,nd' file # delete lines ranging from m to n
# makeblastdb of all 10 protein sets
makeblastdb -in allproteins.fasta -dbtype prot -title "Cyanobacteria 10 taxa genome-wide proteins" -parse_seqids -hash_index -out cyano_prot


# extract sequences from fasta file by ID
extractseq -sequence fasta::"inputfile":"seqid" -outseq "outfile" -auto
# Then run blastp on the sequence
blastp -query "queryfile" -subject "db/fasta" > "outfile"
# Sample evaluation of internal consistency orthogroup prediction
# Shown for inflation value 1.3, orthogroup OG0000519
for protein in {WP_012164843.1,WP_012167703.1,WP_006623986.1,WP_007307728.1,WP_007312090.1,WP_012593649.1,WP_011140679.1,WP_002759324.1,WP_041804180.1,WP_006196060.1,WP_042201566.1,WP_013191945.1,WP_020861411.1}; do
  input="../proteinsets/allproteins.fasta"
  database="../proteinsets/cyano_prot"
  outfile1="seq_${protein}"
  outfile2="iv1.3_OG0000519.txt"
  echo -e "Processing sequence ID $protein"
  extractseq -sequence fasta::$input:$protein -outseq $outfile1 -auto
  blastp -query $outfile1 -db $database >> $outfile2
done