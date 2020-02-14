#!/bin/bash -l

# This script will clone mcclintock into the desired directory and then launch all of the analysis associated with mcclintock.
# Output directory is provided by the user as a command line argument

date=`date +%d_%m_%y`
output_dir=$1"/mcclintock_"$date
mkdir -p $output_dir
dir=`pwd`

git clone https://github.com/bergmanlab/mcclintock.git $output_dir/mcclintock
cd $output_dir/mcclintock
bash install.sh

# Run the reference sequence / TE annotation parameter option analysis on a simulated datasets using sacCer2 reference genome at 3 different sequencing depths
# 4 runs per depth dataset with Repeatmasker TE instances (reference genome; reference genome + consensus TEs; reference genome + reference TE instances; reference genome + consensus TEs + reference TE instances)
# 4 runs per depth dataset with Carr et al TE instances (reference genome; reference genome + consensus TEs; reference genome + reference TE instances; reference genome + consensus TEs + reference TE instances)
cd $dir
cd sacCer2_reference_simulations
echo "bash run_reference_simulations.sh $output_dir" | qsub -V -cwd -l cores=4 -l mem=36 -N mcc_sacCer2_simulations_setup

# Run McClintock on 299 synthetic samples with 1 insertion per sample
cd $dir
cd sacCer2_synthetic_insertion_simulations
echo "bash insertionstest.sh $output_dir" | qsub -V -cwd -l cores=4 -N sacCer2_synthetic_insertion_simulations_setup

# Run McClintock on 299 synthetic samples with 1 insertion per sample on the reverse strand
cd $dir
cd sacCer2_synthetic_insertion_simulations
echo "bash insertionstestrev.sh $output_dir" | qsub -V -cwd -l cores=4 -N sacCer2_synthetic_insertion_simulations_setup

# Run McClintock on the entire McCusker dataset
cd $dir
cd mccusker
echo "bash run_mccusker.sh $output_dir" | qsub -V -cwd -N mccusker_setup
