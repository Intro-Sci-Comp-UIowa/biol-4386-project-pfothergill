#!/bin/bash

run_dir=/Users/pfothergill/
output_dir=/Users/pfothergill/GitHub_Projects/biol-6298-pfothergill/output 
data_dir=/Users/pfothergill/GitHub_Projects/biol-6298-pfothergill/data
#mkdir -p $output_dir/data
mkdir -p $output_dir/qsub_output

# Download the reference genome from UCSC (allows easy browsing of results)
printf "Downloading reference genome...\n\n"
wget -nc http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/chromFa.tar.gz -O $data_dir/chromFa.tar.gz 
tar xvzf $data_dir/chromFa.tar.gz -C $data_dir/
rm $data_dir/chromFa.tar.gz
cat $data_dir/chr*fa $data_dir/2micron.fa > $data_dir/sacCer2.fasta
rm $data_dir/chr*fa $data_dir/2micron.fa

# Download gff locations of reference TE copies
wget -nc http://files.figshare.com/287395/File_S2.txt -O $data_dir/File_S2.txt
awk '{print $3"\treannotate\ttransposable_element\t"$4"\t"$5"\t.\t"$6"\t.\tID="$1;}' $data_dir/File_S2.txt > $data_dir/tmp
sed '1d;$d' $data_dir/tmp > $data_dir/reference_TE_locations.gff
rm $data_dir/File_S2.txt
rm $data_dir/tmp

# The TE families file and consensus TE fasta file are included in this folder > $data_dir/sample_list
# Download list of file locations from EBI
projects=('SRA072302')
for project in "${projects[@]}"
do
	wget -O $data_dir/$project "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$project&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,scientific_name,instrument_model,library_layout,library_source,library_selection,read_count,base_count,experiment_title,fastq_ftp"
	sed 1d $data_dir/$project >> $data_dir/sample_list
done

# Run the pipeline
awk -F'[\t;]' '{print $15}' $data_dir/sample_list > $data_dir/sample_1_urls.txt
#head -5 $data_dir/sample_1_urls.txt > $data_dir/tmp mv $data_dir/tmp $data_dir/sample_1_urls.txt

sample_no=1

while read line
do
	if [ $sample_no -eq 1 ]; then
		qsub -l cores=4 -l long -V -cwd -N mcclintock_mccusker_$sample_no -o $output_dir"/qsub_output/"$sample_no".o" -e $output_dir"/qsub_output/"$sample_no".e" launchmcclintock.sh $line $output_dir $run_dir
		echo -e "$line\tmcclintock$sample_no" > $data_dir/job_key
	else
		qsub -l cores=4 -l long -V -cwd -N mcclintock_mccusker_$sample_no -hold_jid mcclintock_mccusker_"1" -o $output_dir"/qsub_output/"$sample_no".o" -e $output_dir"/qsub_output/"$sample_no".e" launchmcclintock.sh $line $output_dir $run_dir
		echo -e "$line\tmcclintock$sample_no" >> $data_dir/job_key
	fi
	sample_no=$((sample_no+1))
done < $data_dir/sample_1_urls.txt
