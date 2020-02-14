#! /bin/bash -l

samplenumber=$1
run_dir=$2
output_dir=$2"/singleinsertionsims"
sample_name=insertion$samplenumber

# If the sample folder exists this sample must have at least been launched
if [[ -d $output_dir/$sample_name/ ]]
then
	# If the results files exist then nothing needs to be done
	if [[ -f $output_dir/$sample_name/results/$sample_name"_relocate_nonredundant.bed" && -f $output_dir/$sample_name/results/$sample_name"_temp_nonredundant.bed"  && -f $output_dir/$sample_name/results/$sample_name"_ngs_te_mapper_nonredundant.bed" && -f $output_dir/$sample_name/results/$sample_name"_retroseq_nonredundant.bed" && -f $output_dir/$sample_name/results/$sample_name"_telocate_nonredundant.bed" && -f $output_dir/$sample_name/results/$sample_name"_popoolationte_nonredundant.bed" ]]
	then
		echo "$sample_name has already been analysed"
		# If they don't then the run must have been interrupted so files are cleaned up and the run is attempted again
	else
		echo "$sample_name has already been started but must have failed to complete. Cleaning up and retrying"
		rm -r $test_dir/sacCer2/$sample_name

		# Create the genome containing the chromosome with an insertion to sample
		editedchr=`head -1 $output_dir/data/$samplenumber.fasta | cut -d \> -f2`
		editedchr="$editedchr.fa"
		chromosomes=( chrI.fa chrII.fa chrIII.fa chrIV.fa chrV.fa chrVI.fa chrVII.fa chrVIII.fa chrIX.fa chrX.fa chrXI.fa chrXII.fa chrXIII.fa chrXIV.fa chrXV.fa chrXVI.fa chrM.fa 2micron.fa)
		i=0
		for i in {0..17}
		do
			if [[ ${chromosomes[$i]} == $editedchr ]]
			then
				cat $output_dir/data/$samplenumber.fasta >> $output_dir/data/$sample_name"_syntheticsamplegenome.fasta"
			else
				cat $output_dir/data/${chromosomes[$i]} >> $output_dir/data/$sample_name"_syntheticsamplegenome.fasta"
			fi
		done
	
		$run_dir/mcclintock/scripts/fixfastalinelength.pl $output_dir/data/$sample_name"_syntheticsamplegenome.fasta" 80 $output_dir/data/$sample_name"_syntheticsamplegenomefixed.fasta"	
		samtools faidx $output_dir/data/$sample_name"_syntheticsamplegenomefixed.fasta"
		pairs=`awk '{ sum+=$2} END {printf "%.0f", (sum*100)/202}' $output_dir/data/$sample_name"_syntheticsamplegenomefixed.fasta.fai"`
		echo $pairs

		# Create simulated reads for this "sample"
		wgsim -1 101 -2 101 -N $pairs -S 42 -e 0.01 -d 300 -h $output_dir/data/$sample_name"_syntheticsamplegenome.fasta" $output_dir/data/$sample_name"_1.fastq" $output_dir/data/$sample_name"_2.fastq" > $output_dir/data/$sample_name"_wgsimreport"

		# Run the pipeline
		echo "About to launch mcclintock for $sample_name"
		$run_dir/mcclintock/mcclintock.sh -o $output_dir -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $output_dir/data/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 $output_dir/data//$sample_name"_1.fastq" -2 $output_dir/data/$sample_name"_2.fastq" -p 4 -M 35 -i

		# Remove the fastq files and synthetic genome when analysis is complete.
		rm $output_dir/data/$sample_name"_1.fastq" $output_dir/data/$sample_name"_2.fastq" $output_dir/data/$sample_name"_syntheticsamplegenome.fasta"
	fi
else
	# Create the genome containing the chromosome with an insertion to sample
	editedchr=`head -1 $output_dir/data/$samplenumber.fasta | cut -d \> -f2`
	editedchr="$editedchr.fa"
	chromosomes=( chrI.fa chrII.fa chrIII.fa chrIV.fa chrV.fa chrVI.fa chrVII.fa chrVIII.fa chrIX.fa chrX.fa chrXI.fa chrXII.fa chrXIII.fa chrXIV.fa chrXV.fa chrXVI.fa chrM.fa 2micron.fa)
	i=0
	for i in {0..17}
	do
		if [[ ${chromosomes[$i]} == $editedchr ]]
		then
			cat $output_dir/data/$samplenumber.fasta >> $output_dir/data/$sample_name"_syntheticsamplegenome.fasta"
		else
			cat $output_dir/data/${chromosomes[$i]} >> $output_dir/data/$sample_name"_syntheticsamplegenome.fasta"
		fi
	done

	$run_dir/mcclintock/scripts/fixfastalinelength.pl $output_dir/data/$sample_name"_syntheticsamplegenome.fasta" 80 $output_dir/data/$sample_name"_syntheticsamplegenomefixed.fasta"
	samtools faidx $output_dir/data/$sample_name"_syntheticsamplegenomefixed.fasta"  
	pairs=`awk '{ sum+=$2} END {printf "%.0f", (sum*100)/202}' $output_dir/data/$sample_name"_syntheticsamplegenomefixed.fasta.fai"`
	echo $pairs

	# Create simulated reads for this "sample"
	wgsim -1 101 -2 101 -N $pairs -S 42 -e 0.01 -d 300 -h $output_dir/data/$sample_name"_syntheticsamplegenome.fasta" $output_dir/data/$sample_name"_1.fastq" $output_dir/data/$sample_name"_2.fastq" > $output_dir/data/$sample_name"_wgsimreport"

	# Run the pipeline
	echo "About to launch mcclintock for $sample_name"
	$run_dir/mcclintock/mcclintock.sh -o $output_dir -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $output_dir/data/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 $output_dir/data//$sample_name"_1.fastq" -2 $output_dir/data/$sample_name"_2.fastq" -p 4 -M 35 -i

	# Remove the fastq files and synthetic genome when analysis is complete.
	rm $output_dir/data/$sample_name"_1.fastq" $output_dir/data/$sample_name"_2.fastq" $output_dir/data/$sample_name"_syntheticsamplegenome.fasta"
fi
