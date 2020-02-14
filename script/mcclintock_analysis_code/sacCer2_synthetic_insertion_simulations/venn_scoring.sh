# Script to produce matrices for Venn diagrams
# The first argument should be the absolute path to the singleinsertionsims folder
# The second argument should be the absolute path to the output folder location

direction=0
directory=$1
output_directory=$2

> $output_directory"/vennmatrix.tsv"

# Loop once for forward insertions and once for reverse
for run in {1..2}
do
	# Loop through single insertion locations and check for a prediction at this location - if there is one then record the location ID number
	count=1
	while read p
	do
		printf $((count+direction)) >> $output_directory"/vennmatrix.tsv"
		printf "\t" >> $output_directory"/vennmatrix.tsv"
		num=`bedtools window -w 0 -a <(echo "$p") -b <(grep non- $directory"/sacCer2/insertion$count/results/insertion"$count"_ngs_te_mapper_nonredundant.bed") | awk -F'_non' '{print $1}'| awk '{if ($4==$8) print $0}' | wc | awk '{print $1}'`
		if [[ $num > 0 ]]
		then
			printf $((count+direction)) >> $output_directory"/vennmatrix.tsv"
		fi
		printf "\t" >> $output_directory"/vennmatrix.tsv"
		num=`bedtools window -w 0 -a <(echo "$p") -b <(grep non- $directory"/sacCer2/insertion$count/results/insertion"$count"_relocate_nonredundant.bed")  | awk -F'_non' '{print $1}'| awk '{if ($4==$8) print $0}' | wc | awk '{print $1}'`
		if [[ $num > 0 ]]
		then
			printf $((count+direction)) >> $output_directory"/vennmatrix.tsv"
		fi
		printf "\t" >> $output_directory"/vennmatrix.tsv"
		num=`bedtools window -w 0 -a <(echo "$p") -b <(grep non- $directory"/sacCer2/insertion$count/results/insertion"$count"_temp_nonredundant.bed" | grep _sr_)  | awk -F'_non' '{print $1}'| awk '{if ($4==$8) print $0}' | wc | awk '{print $1}'`
		if [[ $num > 0 ]]
		then
			printf $((count+direction)) >> $output_directory"/vennmatrix.tsv"
		fi
		printf "\t" >> $output_directory"/vennmatrix.tsv"
		num=`bedtools window -w 0 -a <(echo "$p") -b <(grep non- $directory"/sacCer2/insertion$count/results/insertion"$count"_temp_nonredundant.bed" | grep _rp_)  | awk -F'_non' '{print $1}'| awk '{if ($4==$8) print $0}' | wc | awk '{print $1}'`
		if [[ $num > 0 ]]
		then
			printf $((count+direction)) >> $output_directory"/vennmatrix.tsv"
		fi
		printf "\t" >> $output_directory"/vennmatrix.tsv"
		num=`bedtools window -w 100 -a <(echo "$p") -b <(grep non- $directory"/sacCer2/insertion$count/results/insertion"$count"_retroseq_nonredundant.bed")  | awk -F'_non' '{print $1}'| awk '{if ($4==$8) print $0}' | wc | awk '{print $1}'`
		if [[ $num > 0 ]]
		then
			printf $((count+direction)) >> $output_directory"/vennmatrix.tsv"
		fi
		printf "\t" >> $output_directory"/vennmatrix.tsv"
		num=`bedtools window -w 0 -a <(echo "$p") -b <(grep non- $directory"/sacCer2/insertion$count/results/insertion"$count"_popoolationte_nonredundant.bed")  | awk -F'_non' '{print $1}'| awk '{if ($4==$8) print $0}' | wc | awk '{print $1}'`
		if [[ $num > 0 ]]
		then
			printf $((count+direction)) >> $output_directory"/vennmatrix.tsv"
		fi
		printf "\t" >> $output_directory"/vennmatrix.tsv"
		num=`bedtools window -w 500 -a <(echo "$p") -b <(grep non- $directory"/sacCer2/insertion$count/results/insertion"$count"_telocate_nonredundant.bed")  | awk -F'_non' '{print $1}'| awk '{if ($4==$8) print $0}' | wc | awk '{print $1}'`
		if [[ $num > 0 ]]
		then
			printf $((count+direction)) >> $output_directory"/vennmatrix.tsv"
		fi
		printf "\n" >> $output_directory"/vennmatrix.tsv"
		count=$((count+1))
	done < $output_directory/singleinsertionlocations.bed
	# Now calculate for the reverse insertions, need to add 299 to count and change directory
	direction=299
	directory=$directory"/../singlerevinsertionsims"
done
