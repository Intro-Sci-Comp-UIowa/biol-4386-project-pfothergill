#! /bin/bash -l

rundir=$1
output_dir=$1"/singlerevinsertionsims"
mkdir -p $output_dir
mkdir -p $output_dir/data
mkdir -p $output_dir/qsub_output

# Download sacCer2 annotation and extract tRNA coordinates
wget -q http://downloads.yeastgenome.org/curation/chromosomal_feature/archive/saccharomyces_cerevisiae_R61-1-1_20080607.gff.gz -O $output_dir/data/saccharomyces_cerevisiae_R61-1-1_20080607.gff.gz
gunzip -c $output_dir/data/saccharomyces_cerevisiae_R61-1-1_20080607.gff.gz > $output_dir/data/saccharomyces_cerevisiae_R61-1-1_20080607.gff
rm $output_dir/data/saccharomyces_cerevisiae_R61-1-1_20080607.gff.gz
awk '{if ($3=="tRNA") print $0}' $output_dir/data/saccharomyces_cerevisiae_R61-1-1_20080607.gff > $output_dir/data/tRNA_coordinates.gff
rm $output_dir/data/saccharomyces_cerevisiae_R61-1-1_20080607.gff

# Download gff locations of reference TE copies
wget -nc -q http://files.figshare.com/287395/File_S2.txt -O $output_dir/data/File_S2.txt
awk '{print $3"\treannotate\ttransposable_element\t"$4"\t"$5"\t.\t"$6"\t.\tID="$1;}' $output_dir/data/File_S2.txt > $output_dir/data/tmp
sed '1d;$d' $output_dir/data/tmp > $output_dir/data/reference_TE_locations.gff
rm $output_dir/data/File_S2.txt
rm $output_dir/data/tmp

# Download and extract sacCer2 chromosomes
wget -nc -q http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/chromFa.tar.gz -O $output_dir/data/chromFa.tar.gz
tar xvzf $output_dir/data/chromFa.tar.gz -C $output_dir/data/
rm $output_dir/data/chromFa.tar.gz
cat $output_dir/data/chr*fa $output_dir/data/2micron.fa > $output_dir/data/sacCer2.fasta

# Index reference to get chr lengths
samtools faidx $output_dir/data/sacCer2.fasta

# Edit chr name for mito to match UCSC
sed -i 's/chrMito/chrM/g' $output_dir/data/tRNA_coordinates.gff

# Get start and end chromosome sections. Every third set of coordinates will be a TY3 and so uses a different insertion location.
awk 'FNR==NR{a[$1]=$2;next}{if ($7=="+") {if ((FNR+1)%4==0) {print $1"\t0\t"$4-13"\t"$1"\n"$1"\t"$4-18"\t"a[$1]"\t"$1} else {print $1"\t0\t"$4-196"\t"$1"\n"$1"\t"$4-201"\t"a[$1]"\t"$1}} else {if ((FNR+1)%4==0) {print $1"\t0\t"$4+17"\t"$1"\n"$1"\t"$4+12"\t"a[$1]"\t"$1} else {print $1"\t0\t"$4+200"\t"$1"\n"$1"\t"$4+195"\t"a[$1]"\t"$1}}}' $output_dir/data/sacCer2.fasta.fai $output_dir/data/tRNA_coordinates.gff > $output_dir/data/insertions.bed

# Extract start and end fasta sequence now with a TSD
bedtools getfasta -name -fi $output_dir/data/sacCer2.fasta -bed $output_dir/data/insertions.bed -fo $output_dir/data/fragments.fasta

TElines=("" "80,199" "508,627" "400,507" "200,290")

# Insert a TE into the chromosome
TE=1
j=0

for i in {1..299}
do
	start=$((i*2+$j))
	id=$((start-1))
	j=$((j+2))
	end=$((i*2+$j))

	line=$((i*2-1))
	awk -v TEnum="$TE" -v line="$line" '{if(NR==line) {end=$3; getline; print $1"\t"$2"\t"end"\tTY"TEnum} }' $output_dir/data/insertions.bed > $output_dir/data/$i".bed"

	sed -n "${TElines[$TE]} p" $rundir/mcclintock/test/sac_cer_TE_seqs.fasta > $rundir/mcclintock/test/tempTE.fasta
	$rundir/mcclintock/scripts/fixfastalinelength.pl $rundir/mcclintock/test/tempTE.fasta 80 $rundir/mcclintock/test/fixedtempTE.fasta	
	revFasta $rundir/mcclintock/test/fixedtempTE.fasta > $rundir/mcclintock/test/revTE.fasta
	
	cat <(sed -n "$id,$start p" $output_dir/data/fragments.fasta) <(tail -n +2 $rundir/mcclintock/test/revTE.fasta)  <(sed -n "$end p" $output_dir/data/fragments.fasta) > $output_dir/data/$i".fasta"

	# Loop through the 4 different mobile TEs in sac_cer_TE_seqs.fasta
	TE=$((TE+1))
	if [[ $TE -eq 5 ]]
	then
		TE=1
	fi
done

# Launch the first sample
echo "bash launchmcclintockrev.sh 1 $rundir" | qsub -V -cwd -l cores=4 -l long -N insertion_1 -o $output_dir"/qsub_output/insertion_1.o" -e $output_dir"/qsub_output/insertion_1.e"

for i in {2..299}
do
	echo "bash launchmcclintockrev.sh $i $rundir" | qsub -V -cwd -l cores=4 -l long -N insertion_$i -hold_jid insertion_1 -o $output_dir"/qsub_output/insertion_"$i".o" -e $output_dir"/qsub_output/insertion_"$i".e"
done

