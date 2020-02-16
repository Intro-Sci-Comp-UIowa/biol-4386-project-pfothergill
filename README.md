# Reference(s)
  
- **Nelson, M. G., R. S. Linheiro, and C. M. Bergman, 2017  McClintock: An Integrated Pipeline for Detecting Transposable Element Insertions in Whole-Genome Shotgun Sequencing Data. G3 (Bethesda) 7: 2763-2778**

- Strope P. K., Skelly D. A., Kozmin S. G., Mahadevan G., Stone E. A., et al., 2015  The 100-genomes strains, an S. cerevisiae resource that illuminates its natural phenotypic and genotypic variation and emergence as an opportunistic pathogen. Genome Res. 25: 762-774.  

# Introduction 
  
- Using the [mcclintock pipeline](https://github.com/bergmanlab/mcclintock) I am going to recreate [Figure 5](https://www.g3journal.org/content/ggg/7/8/2763/F5.large.jpg?width=800&height=600&carousel=1) from the [paper](https://www.g3journal.org/content/7/8/2763) that's used to help setup and test the pipeline with real world data from S. cerevisiae. The genomes are to be downloaded from the [Strope et al](https://genome.cshlp.org/content/25/5/762?ijkey=704a2c3eaf47b3364b45fabb99243292adddcd05&keytype2=tf_ipsecsha) paper. This pipeline finds transposable element insertions/deletions in reference and non-reference whole-genome shotgun sequences.
 
# Figure  
  
[<img src="https://www.g3journal.org/content/ggg/7/8/2763/F5.large.jpg?width=800&height=600&carousel=1">](https://www.g3journal.org/content/ggg/7/8/2763/F5.large.jpg?width=800&height=600&carousel=1)  

# Outline of approaches
1. Setup Conda Environment (For use of python3 tools during pipeline run)
	- Follow the first and second installation instructions (1 & 2) on GitHub README
		- Download Miniconda (python3)
		- Update Conda
2. Setup McClintock Pipeline
	- Clone McClintock GitHub 
		- Using SSH: 
			- $ ```git clone git@github.com:bergmanlab/mcclintock.git```
		- Using HTTPS:
			- $ ```git clone https://github.com/bergmanlab/mcclintock.git```
	- Follow the third and fourth installation instruction (3 & 4) on GitHub README
		- Install Conda dependencies by running provided env_installer.sh script
		- Install the pipeline using the provided install.sh script
3. Download all my fastq files (93 strains x 2 reads/strain)
	- I need about 15.5 days to download all materials as each fastq file takes about 2 hours to fully download and I have 186 of them
	- (OPTIONAL) Use a script to install all my fastq files (93 x 2) in local McClintock folder in Argon
4. Write a script to make 93 consecutive runs, one for each strain, that outputs into a specified "output" folder
5. Plotting
	- (OPTIONAL) Use SeqPlots to make a boxplot of my data by running the program and telling the program the right place to find all my data
	- If above doesn't work correctly, manually plot all 558 data points corresponding to each program ran on each stain. Tedious but possibly my only option
		- Count lines for each BED file of each strain and map that on the boxplot for the corresponding program. Do this for all 558 BED files
