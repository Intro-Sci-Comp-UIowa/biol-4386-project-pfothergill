# Reference
  
- Nelson, M. G., R. S. Linheiro, and C. M. Bergman, 2017 McClintock: An Integrated Pipeline for Detecting Transposable Element Insertions in Whole-Genome Shotgun Sequencing Data. G3 (Bethesda) 7: 2763–2778.  
  
# Introduction  
  
- Using the [mcclintock pipeline](https://github.com/bergmanlab/mcclintock) I am going to recreate [Figure 5](https://www.g3journal.org/content/ggg/7/8/2763/F5.large.jpg?width=800&height=600&carousel=1) from the [paper](https://www.g3journal.org/content/7/8/2763) that's used to help setup and test the pipeline. This pipeline finds transposable element insertions/deletions in reference and non-reference whole-genome shotgun sequences.
 
# Figure  
  
[<img src="https://www.g3journal.org/content/ggg/7/8/2763/F5.large.jpg?width=800&height=600&carousel=1">](https://www.g3journal.org/content/ggg/7/8/2763/F5.large.jpg?width=800&height=600&carousel=1)  

# Outline of approaches

1. Install miniconda
2. Add conda channels
3. Clone bergmanlab/mcclintock git repository
4. cd into the repo and run the environment install and pipeline installer
5. Run the test dataset through the pipeline using bash/terminal
6. Find the output directory and use a program to make a box plot of the results in the directory.