#/bin/bash

# Search for the word "non-reference in each output strain for the corresponding program"
grep -c "non-reference" ../output/SRA072302/SRR800*/results/*ngs_te_mapper*.bed > ngs_te_mapper_full_data.txt
grep -c "non-reference" ../output/SRA072302/SRR800*/results/*popoolationte*.bed > popoolationte_full_data.txt
grep -c "non-reference" ../output/SRA072302/SRR800*/results/*relocate*.bed > relocate_full_data.txt
grep -c "non-reference" ../output/SRA072302/SRR800*/results/*retroseq*.bed > retroseq_full_data.txt
grep -c "non-reference" ../output/SRA072302/SRR800*/results/*telocate*.bed > telocate_full_data.txt
grep -c "non-reference" ../output/SRA072302/SRR800*/results/*temp*.bed > temp_full_data.txt

# Formating stuff, we want to get the second field after the ":" from the output of the commands above
cut -d ":" -f 2 ngs_te_mapper_full_data.txt | sort -n > ngs_te_mapper_sorted_freq.txt
cut -d ":" -f 2 popoolationte_full_data.txt | sort -n > popoolationte_sorted_freq.txt
cut -d ":" -f 2 relocate_full_data.txt | sort -n > relocate_sorted_freq.txt
cut -d ":" -f 2 retroseq_full_data.txt | sort -n > retroseq_sorted_freq.txt
cut -d ":" -f 2 telocate_full_data.txt | sort -n > telocate_sorted_freq.txt
cut -d ":" -f 2 temp_full_data.txt | sort -n > temp_sorted_freq.txt

# Next create basename that will later reside at the top of the csv file used in R
for f in *_sorted_freq.txt; do echo $(basename $f _sorted_freq.txt); done > program_basenames.txt

# formatting for the csv file. All intermediate files will be removed at the end and you'll be left
# with the combined.txt (even though it's a csv, it has the extension txt and will work!)
awk 'ORS=NR?" ":"\n"' program_basenames.txt | awk ' { print $1, $3, $6, $4, $2, $5  } ' > horizontal_basenames.txt
paste ngs_te_mapper_sorted_freq.txt \
relocate_sorted_freq.txt \
temp_sorted_freq.txt \
retroseq_sorted_freq.txt \
popoolationte_sorted_freq.txt \
telocate_sorted_freq.txt > all_sorted.txt
cat horizontal_basenames.txt all_sorted.txt | column -t > tabular_combined.txt
awk ' { print $1, $2, $3, $4, $5, $6 } ' tabular_combined.txt | sed 's/ /,/g' > combined.txt
rm all_sorted.txt \
tabular_combined.txt \
program_basenames.txt \
ngs* \
horizontal_basenames.txt \
popoolationte* \
relocate* \
retroseq* \
telocate* \
temp* \

# run the R script that will make the recreated box-plot
chmod +x make_plot.R
R --no-save < make_plot.R
mkdir -p ../output/recreated_image
mv recreated_boxplot.png ../output/recreated_image