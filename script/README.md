# Instructions for every script
1. Getting the mcclintock_analysis_code directory - useful because it contains scripts for running the entire experiment (from obtaining fastq files, reference genomes, gff file of TE locations, running the pipeline, etc.)
 	```
	$ cd <project directory root>/script
	$ wget https://www.g3journal.org/highwire/filestream/485221/field_highwire_adjunct_files/2/FileS3.zip
	$ unzip FileS3.zip
	```

	- In the resulting files, I had to make changes to both the mccusker/run_mccusker.sh and mccusker/launchmcclintock.sh. The changes (for both files) is changing every single "$output_dir/data" -> "$data_dir" and at the top of both files, setting "data_dir=<project directory root>/data". I also changed run_dir in the mccusker/run_mccusker.sh file to the absolute path to my ***INSTALLED*** McClintock pipeline. In my case, the mcclintock pipeline was installed to to my $HOME directory so I changed run_dir to be run_dir=/Users/pfothergill instead of $1. I also changed the output directory to output_dir=<project directory root>/output & commented out "mkdir -p $output_dir/data"
