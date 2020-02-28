# Instructions for data reproduction

- I have put a few scripts in the ../script/mcclintock_analysis_code/mccusker titled:
	1. run_mccusker.sh
	2. launchmcclintock.sh
- The first one will initialize information such as the run directory, output directory, where to wget the strain fastq's, etc.
- At the end, run_mccusker.sh will call launchmcclintock.sh which then will use the information from run_mccusker.sh to run through the mccclintock pipeline.
- I have started the run_mccusker.sh script, however, if there isn't anything in this directory other than this README file, it means it is still running and I haven't been able to update it with all the fastq's from all of 93 strains and both of a reference genome (will be titled sacCer2.fasta) and gff locations of reference TE copies (will be titled reference_TE_locations.gff).
- Data will be coming from:
	1. http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$project&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,scientific_name,instrument_model,library_layout,library_source,library_selection,read_count,base_count,experiment_title,fastq_ftp
		- Accounts for genomes of all 93 strains
	2. http://files.figshare.com/287395/File_S2.txt
		- GFF locations of reference TE copies
	3. http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/chromFa.tar.gz
		- Contains 2 files - once I untar the file, the 2 files will be concatenated together

# To obtain all data and run the pipeline:
```
$ cd <project root directory>/script/mcclintock_analysis_code/mccusker
$ ./run_mccusker
```
- This will work, even for HPC (as long as qsub is used for your HPC) because the last part of the run_mccusker.sh script will make qsub calls with all the newly formed data (from the first part of the script -- wget commands)

# To obtain fastq files only (no pipeline run):
```
$ cd data
$ ./obtain_fastq_files.sh 
```

