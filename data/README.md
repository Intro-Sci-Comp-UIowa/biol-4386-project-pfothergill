# Instructions for data reproduction

- I have put a few scripts in the ../script/mcclintock_analysis_code/mccusker titled:
	1. run_mccusker.sh
	2. launchmcclintock.sh
- The first one will initialize information such as the run directory, output directory, where to wget the strain fastq's, etc.
- At the end, run_mccusker.sh will call launchmcclintock.sh which then will use the information from run_mccusker.sh to run through the mccclintock pipeline.
- I have started the run_mccusker.sh script, however, if there isn't anything in this directory other than this README file, it means it is still running and I haven't been able to update it with all the fastq's from all of 93 strains and both of a reference genome (will be titled sacCer2.fasta) and gff locations of reference TE copies (will be titled reference_TE_locations.gff).
