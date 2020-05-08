# Reference(s)
  
- **Nelson, M. G., R. S. Linheiro, and C. M. Bergman, 2017  McClintock: An Integrated Pipeline for Detecting Transposable Element Insertions in Whole-Genome Shotgun Sequencing Data. G3 (Bethesda) 7: 2763-2778**

- Strope P. K., Skelly D. A., Kozmin S. G., Mahadevan G., Stone E. A., et al., 2015  The 100-genomes strains, an S. cerevisiae resource that illuminates its natural phenotypic and genotypic variation and emergence as an opportunistic pathogen. Genome Res. 25: 762-774.  

# Introduction 
  
- Using the [mcclintock pipeline](https://github.com/bergmanlab/mcclintock) I am going to recreate [Figure 5](https://www.g3journal.org/content/ggg/7/8/2763/F5.large.jpg?width=800&height=600&carousel=1) from the [paper](https://www.g3journal.org/content/7/8/2763) that's used to help setup and test the pipeline with real world data from S. cerevisiae. The genomes are to be downloaded from the [Strope et al](https://genome.cshlp.org/content/25/5/762?ijkey=704a2c3eaf47b3364b45fabb99243292adddcd05&keytype2=tf_ipsecsha) paper. This pipeline finds transposable element insertions/deletions in reference and non-reference whole-genome shotgun sequences.
 
# Figure  
  
[<img src="https://www.g3journal.org/content/ggg/7/8/2763/F5.large.jpg?width=800&height=600&carousel=1">](https://www.g3journal.org/content/ggg/7/8/2763/F5.large.jpg?width=800&height=600&carousel=1)  

# How to run (Installation steps 1-4 are taken directly from the original [bergmanlab](https://github.com/bergmanlab/mcclintock) GitHub page. See that page for more details about installation, testing and each program used in the pipeline)

## Installation

1. Install python3 miniconda (miniconda is a lightweight installer for the conda package manager).
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $HOME//miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda # silent mode
echo "export PATH=\$PATH:\$HOME/miniconda/bin" >> $HOME/.bashrc # add to .bashrc
source $HOME/.bashrc
```

2. Update conda and set up your conda channels (this lets conda know where to find packages)
```
conda update conda
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

3. Install dependencies with conda by clone the McClintock repository, cd into the project directory and run the script install.sh with no arguments, the installation might take a while to finish.
```
git clone git@github.com:bergmanlab/mcclintock.git
cd mcclintock
sh env_install.sh
```

4. Install McClintock by run the script install.sh with no arguments.
```
sh install.sh
```

install.sh will download and unpack all of the TE detection pipelines and check that the required dependencies are available in your path. Missing dependencies will be reported and you must install or make sure these are available to run the full pipeline.

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
3. (OPTIONAL) Download all fastq files (93 strains x 2 reads/strain) individually or used as a scipt
	- This is done automatically within script/mcclintock_analysis_code/run_mccusker.sh so only do this if you want only the reads used in this study
	- Plan accordingly, individually, they could take upwards of 2 weeks to obtain as each fastq file takes about ~2 hours to fully download and there are 186 of them
5. Plotting
	- (OPTIONAL) Use SeqPlots to make a boxplot of my data by running the program and telling the program the right place to find all my data
	- If above doesn't work correctly, manually plot all 558 data points corresponding to each program ran on each stain. Tedious but possibly the only option. See the analysis/ directory for more information
		- Count lines for each BED file of each strain and map on the boxplot in R for the corresponding program. Do this for all 558 BED files
