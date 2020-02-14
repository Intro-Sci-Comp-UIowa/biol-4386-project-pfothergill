#!/bin/bash -l

run_dir=$1
output_dir=$1"/sacCer2_reference_simulations"
mkdir $output_dir
mkdir $output_dir/data
mkdir $output_dir/qsub_output

# Download the reference genome from UCSC (allows easy browsing of results)
printf "Downloading reference genome...\n\n"
wget -nc -q http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/chromFa.tar.gz -O $output_dir/data/chromFa.tar.gz
tar xvzf $output_dir/data/chromFa.tar.gz -C $output_dir/data/
rm $output_dir/data/chromFa.tar.gz
# Combine the chromosomes together
cat $output_dir/data/chr*fa $output_dir/data/2micron.fa > $output_dir/data/sacCer2.fasta
rm $output_dir/data/chr*fa $output_dir/data/2micron.fa

# Download gff locations of reference TE copies
wget -nc -q http://files.figshare.com/287395/File_S2.txt -O $output_dir/data/File_S2.txt
awk '{print $3"\treannotate\ttransposable_element\t"$4"\t"$5"\t.\t"$6"\t.\tID="$1;}' $output_dir/data/File_S2.txt > $output_dir/data/tmp
sed '1d;$d' $output_dir/data/tmp > $output_dir/data/reference_TE_locations.gff
rm $output_dir/data/File_S2.txt
rm $output_dir/data/tmp

printf "\nCreate Simulated fastq files...\n\n"
for i in `seq 1 100`
do
	wgsim -1 101 -2 101 -d 300 -N 6021285 -S $i -e 0.01 -h $output_dir/data/sacCer2.fasta $output_dir"/data/100X"$i"simulation_1.fastq" $output_dir"/data/100X"$i"simulation_2.fastq" > $output_dir"/data/wgsim_mutation_report_"$i".txt"
done

# Subsample using seqtk to produce lower coverage datasets.
for i in `seq 1 100`
do
	seqtk seq -s $i -f 0.1 $output_dir"/data/100X"$i"simulation_1.fastq" > $output_dir"/data/10X"$i"simulation_1.fastq"
	seqtk seq -s $i -f 0.1 $output_dir"/data/100X"$i"simulation_2.fastq" > $output_dir"/data/10X"$i"simulation_2.fastq"
done

# The TE families file and consensus TE fasta file are included with mcclintock

# Run the pipeline for all combinations of reference genome and simulated depths of coverage
cd ..

echo "bash $run_dir/mcclintock/mcclintock.sh -o $output_dir/normalreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $output_dir/data/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 $output_dir/data/10X1simulation_1.fastq -2 $output_dir/data/10X1simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N normalref101 -o $output_dir/qsub_output/10X1normalref.o -e $output_dir/qsub_output/10X1normalref.e
echo "bash  $run_dir/mcclintock/mcclintock.sh -C -o $output_dir/consensusreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $output_dir/data/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 $output_dir/data/10X1simulation_1.fastq -2 $output_dir/data/10X1simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N conref101 -o $output_dir/qsub_output/10X1consensusref.o -e $output_dir/qsub_output/10X1consensusref.e
echo "bash  $run_dir/mcclintock/mcclintock.sh -R -o $output_dir/refcopiesreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $output_dir/data/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 $output_dir/data/10X1simulation_1.fastq -2 $output_dir/data/10X1simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N refref101 -o $output_dir/qsub_output/10X1referenceref.o -e $output_dir/qsub_output/10X1referenceref.e
echo "bash  $run_dir/mcclintock/mcclintock.sh -C -R -o $output_dir/fullreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $output_dir/data/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 $output_dir/data/10X1simulation_1.fastq -2 $output_dir/data/10X1simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N fullref101 -o $output_dir/qsub_output/10X1fullref.o -e $output_dir/qsub_output/10X1fullref.e

echo "bash  $run_dir/mcclintock/mcclintock.sh -o $output_dir/RMnormalreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -1 $output_dir/data/10X1simulation_1.fastq -2 $output_dir/data/10X1simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N RMnormalref101 -o $output_dir/qsub_output/10X1RMnormalref.o -e $output_dir/qsub_output/10X1RMnormalref.e
echo "bash  $run_dir/mcclintock/mcclintock.sh -C -o $output_dir/RMconsensusreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -1 $output_dir/data/10X1simulation_1.fastq -2 $output_dir/data/10X1simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N RMconref101 -o $output_dir/qsub_output/10X1RMconsensusref.o -e $output_dir/qsub_output/10X1RMconsensusref.e
echo "bash  $run_dir/mcclintock/mcclintock.sh -R -o $output_dir/RMrefcopiesreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -1 $output_dir/data/10X1simulation_1.fastq -2 $output_dir/data/10X1simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N RMrefref101 -o $output_dir/qsub_output/10X1RMreferenceref.o -e $output_dir/qsub_output/10X1RMreferenceref.e
echo "bash  $run_dir/mcclintock/mcclintock.sh -C -R -o $output_dir/RMfullreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -1 $output_dir/data/10X1simulation_1.fastq -2 $output_dir/data/10X1simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N RMfullref101 -o $output_dir/qsub_output/10X1RMfullref.o -e $output_dir/qsub_output/10X1RMfullref.e



for i in `seq 2 100`
do
	echo "bash $run_dir/mcclintock/mcclintock.sh -o $output_dir/normalreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $output_dir/data/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 "$output_dir"/data/10X"$i"simulation_1.fastq -2 "$output_dir"/data/10X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N normalref10$i -hold_jid normalref101 -o $output_dir"/qsub_output/10X"$i"normalref.o" -e $output_dir"/qsub_output/10X"$i"normalref.e"
	echo "bash  $run_dir/mcclintock/mcclintock.sh -C -o $output_dir/consensusreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $output_dir/data/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 "$output_dir"/data/10X"$i"simulation_1.fastq -2 "$output_dir"/data/10X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N conref10$i -hold_jid conref101 -o $output_dir"/qsub_output/10X"$i"consensusref.o" -e $output_dir"/qsub_output/10X"$i"consensusref.e"
	echo "bash  $run_dir/mcclintock/mcclintock.sh -R -o $output_dir/refcopiesreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $output_dir/data/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 "$output_dir"/data/10X"$i"simulation_1.fastq -2 "$output_dir"/data/10X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N refref10$i -hold_jid refref101 -o $output_dir"/qsub_output/10X"$i"referenceref.o" -e $output_dir"/qsub_output/10X"$i"referenceref.e"
	echo "bash  $run_dir/mcclintock/mcclintock.sh -C -R -o $output_dir/fullreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $output_dir/data/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 "$output_dir"/data/10X"$i"simulation_1.fastq -2 "$output_dir"/data/10X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N fullref10$i -hold_jid fullref101 -o $output_dir"/qsub_output/10X"$i"fullref.o" -e $output_dir"/qsub_output/10X"$i"fullref.e"

	echo "bash  $run_dir/mcclintock/mcclintock.sh -o $output_dir/RMnormalreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -1 "$output_dir"/data/10X"$i"simulation_1.fastq -2 "$output_dir"/data/10X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N RMnormalref10$i -hold_jid RMnormalref101 -o $output_dir"/qsub_output/10X"$i"RMnormalref.o" -e $output_dir"/qsub_output/10X"$i"RMnormalref.e"
	echo "bash  $run_dir/mcclintock/mcclintock.sh -C -o $output_dir/RMconsensusreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -1 "$output_dir"/data/10X"$i"simulation_1.fastq -2 "$output_dir"/data/10X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N RMconref10$i -hold_jid RMconref101 -o $output_dir"/qsub_output/10X"$i"RMconsensusref.o" -e $output_dir"/qsub_output/10X"$i"RMconsensusref.e"
	echo "bash  $run_dir/mcclintock/mcclintock.sh -R -o $output_dir/RMrefcopiesreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -1 "$output_dir"/data/10X"$i"simulation_1.fastq -2 "$output_dir"/data/10X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N RMrefref10$i -hold_jid RMrefref101 -o $output_dir"/qsub_output/10X"$i"RMreferenceref.o" -e $output_dir"/qsub_output/10X"$i"RMreferenceref.e"
	echo "bash  $run_dir/mcclintock/mcclintock.sh -C -R -o $output_dir/RMfullreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -1 "$output_dir"/data/10X"$i"simulation_1.fastq -2 "$output_dir"/data/10X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N RMfullref10$i -hold_jid fullref101 -o $output_dir"/qsub_output/10X"$i"RMfullref.o" -e $output_dir"/qsub_output/10X"$i"RMfullref.e"
done

for i in `seq 1 100`
do
	echo "bash $run_dir/mcclintock/mcclintock.sh -o $output_dir/normalreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $output_dir/data/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 "$output_dir"/data/100X"$i"simulation_1.fastq -2 "$output_dir"/data/100X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N normalref100$i -hold_jid normalref101 -o $output_dir"/qsub_output/100X"$i"normalref.o" -e $output_dir"/qsub_output/100X"$i"normalref.e"
	echo "bash  $run_dir/mcclintock/mcclintock.sh -C -o $output_dir/consensusreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $output_dir/data/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 "$output_dir"/data/100X"$i"simulation_1.fastq -2 "$output_dir"/data/100X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N conref100$i -hold_jid conref101 -o $output_dir"/qsub_output/100X"$i"consensusref.o" -e $output_dir"/qsub_output/100X"$i"consensusref.e"
	echo "bash  $run_dir/mcclintock/mcclintock.sh -R -o $output_dir/refcopiesreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $output_dir/data/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 "$output_dir"/data/100X"$i"simulation_1.fastq -2 "$output_dir"/data/100X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N refref100$i -hold_jid refref101 -o $output_dir"/qsub_output/100X"$i"referenceref.o" -e $output_dir"/qsub_output/100X"$i"referenceref.e"
	echo "bash  $run_dir/mcclintock/mcclintock.sh -C -R -o $output_dir/fullreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -g $output_dir/data/reference_TE_locations.gff -t $run_dir/mcclintock/test/sac_cer_te_families.tsv -1 "$output_dir"/data/100X"$i"simulation_1.fastq -2 "$output_dir"/data/100X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N fullref100$i -hold_jid fullref101 -o $output_dir"/qsub_output/100X"$i"fullref.o" -e $output_dir"/qsub_output/100X"$i"fullref.e"

	echo "bash  $run_dir/mcclintock/mcclintock.sh -o $output_dir/RMnormalreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -1 "$output_dir"/data/100X"$i"simulation_1.fastq -2 "$output_dir"/data/100X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N RMnormalref100$i -hold_jid RMnormalref101 -o $output_dir"/qsub_output/100X"$i"RMnormalref.o" -e $output_dir"/qsub_output/100X"$i"RMnormalref.e"
	echo "bash  $run_dir/mcclintock/mcclintock.sh -C -o $output_dir/RMconsensusreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -1 "$output_dir"/data/100X"$i"simulation_1.fastq -2 "$output_dir"/data/100X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80 -N RMconref100$i -hold_jid RMconref101 -o $output_dir"/qsub_output/100X"$i"RMconsensusref.o" -e $output_dir"/qsub_output/100X"$i"RMconsensusref.e"
	echo "bash  $run_dir/mcclintock/mcclintock.sh -R -o $output_dir/RMrefcopiesreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -1 "$output_dir"/data/100X"$i"simulation_1.fastq -2 "$output_dir"/data/100X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80  -N RMrefref100$i -hold_jid RMrefref101 -o $output_dir"/qsub_output/100X"$i"RMreferenceref.o" -e $output_dir"/qsub_output/100X"$i"RMreferenceref.e"
	echo "bash  $run_dir/mcclintock/mcclintock.sh -C -R -o $output_dir/RMfullreference -r $output_dir/data/sacCer2.fasta -c $run_dir/mcclintock/test/sac_cer_TE_seqs.fasta -1 "$output_dir"/data/100X"$i"simulation_1.fastq -2 "$output_dir"/data/100X"$i"simulation_2.fastq -p 4 -M 75 -i" | qsub -cwd -V -l cores=4  -l mem=80  -N RMfullref100$i -hold_jid RMfullref101 -o $output_dir"/qsub_output/100X"$i"RMfullref.o" -e $output_dir"/qsub_output/100X"$i"RMfullref.e"

done

