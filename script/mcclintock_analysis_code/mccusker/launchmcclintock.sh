#!/bin/bash

urltmp=${1%%_1.f*}
url=${urltmp%%_2.f*}
sample_name=${url##*/}
output_dir=$2
run_dir=$3 
data_dir=/Users/pfothergill/GitHub_Projects/biol-6298-pfothergill/data

# If the sample folder exists this sample must have at least been launched
if [ -d $output_dir/$sample_name/ ]; then
	# If the results files exist then nothing needs to be done
	if [[ -f $output_dir/$sample_name/results/$sample_name"_relocate_nonredundant.bed" && -f $output_dir/$sample_name/results/$sample_name"_temp_nonredundant.bed" && -f $output_dir/$sample_name/results/$sample_name"_ngs_te_mapper_nonredundant.bed" && -f $output_dir/$sample_name/results/$sample_name"_retroseq_nonredundant.bed" && -f $output_dir/$sample_name/results/$sample_name"_telocate_nonredundant.bed" && -f $output_dir/$sample_name/results/$sample_name"_popoolationte_nonredundant.bed" ]]
	then
		echo "$sample_name has already been analysed"
		# If they don't then the run must have been interrupted so files are cleaned up and the run is 
		# attempted again
	else
		echo "$sample_name has already been started but must have failed to complete. Cleaning up and retrying"
		rm -r $test_dir/sacCer2/$sample_name

		# wget the fastq files
		while ! wget -c -q $url"_1.fastq.gz" -O $data_dir/$sample_name"_1.fastq.gz"; do sleep 5; done
		while ! wget -c -q $url"_2.fastq.gz" -O $data_dir/$sample_name"_2.fastq.gz"; do sleep 5; done
		
		# Uncompress the fastq file(s)
		gunzip $data_dir/$sample_name"_1.fastq.gz"
		gunzip $data_dir/$sample_name"_2.fastq.gz"
		
		# Lauch the McClintock Pipeline for the sample that was just downloaded
		echo "About to launch mcclintock for $sample_name"
		bash $run_dir/mcclintock/mcclintock.sh -o $output_dir -r $data_dir/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $data_dir/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 $data_dir/$sample_name"_1.fastq" -2 $data_dir/$sample_name"_2.fastq" -i -p 4
	
		# Remove the fastq files when analysis is complete. rm $data_dir/$sample_name"_1.fastq" 
		# $data_dir/$sample_name"_2.fastq"
	fi 
else
	# wget the fastq files
	while ! wget -c -q $url"_1.fastq.gz" -O $data_dir/$sample_name"_1.fastq.gz"; do sleep 5; done
	while ! wget -c -q $url"_2.fastq.gz" -O $data_dir/$sample_name"_2.fastq.gz"; do sleep 5; done
	
	# Uncompress the fastq file(s)
	gunzip $data_dir/$sample_name"_1.fastq.gz"
	gunzip $data_dir/$sample_name"_2.fastq.gz"
		
	echo "About to launch mcclintock for $sample_name"
	bash $run_dir/mcclintock/mcclintock.sh -o $output_dir -r $data_dir/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $data_dir/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 $data_dir/$sample_name"_1.fastq" -2 $data_dir/$sample_name"_2.fastq" -i -p 4
	
	# Remove the fastq files when analysis is complete. rm $data_dir/$sample_name"_1.fastq" 
	# $data_dir/$sample_name"_2.fastq"
fi
