# To **analyze** the data, do the following!

```
$ ./run_analysis.sh
```

## This will do the following (you can run each command seperately if you do not want to run a script for the whole process!):
### First, change into the analysis directory (you are probably here already if you are reading this!)
```
$ cd analysis
```

### Go through and make a file from the output directory that contains the amount of insertions. Insertions can be counted because each programs file in every strain contains exactly 1 insertion per line of the file. So, count the lines of the file (only count "non-reference" lines) and you have the non-reference insertion total. Count all insertions for each specific file of every strain and clump them together by program type.
```
$ grep -c "non-reference" ../output/SRA072302/SRR800*/results/*ngs_te_mapper*.bed > ngs_te_mapper_full_data.txt
$ grep -c "non-reference" ../output/SRA072302/SRR800*/results/*popoolationte*.bed > popoolationte_full_data.txt
$ grep -c "non-reference" ../output/SRA072302/SRR800*/results/*relocate*.bed > relocate_full_data.txt
$ grep -c "non-reference" ../output/SRA072302/SRR800*/results/*retroseq*.bed > retroseq_full_data.txt
$ grep -c "non-reference" ../output/SRA072302/SRR800*/results/*telocate*.bed > telocate_full_data.txt
$ grep -c "non-reference" ../output/SRA072302/SRR800*/results/*temp*.bed > temp_full_data.txt
```

### Sort each new file and only include the amount of insertion per strain per program (insertions/strain/program)
```
$ cut -d ":" -f 2 ngs_te_mapper_full_data.txt | sort -n > ngs_te_mapper_sorted_freq.txt
$ cut -d ":" -f 2 popoolationte_full_data.txt | sort -n > popoolationte_sorted_freq.txt
$ cut -d ":" -f 2 relocate_full_data.txt | sort -n > relocate_sorted_freq.txt
$ cut -d ":" -f 2 retroseq_full_data.txt | sort -n > retroseq_sorted_freq.txt
$ cut -d ":" -f 2 telocate_full_data.txt | sort -n > telocate_sorted_freq.txt
$ cut -d ":" -f 2 temp_full_data.txt | sort -n > temp_sorted_freq.txt
```

### Make a text file with all the program's basename on individual lines
```
$ for f in *_sorted_freq.txt; do echo $(basename $f _sorted_freq.txt); done > program_basenames.txt
```

### To get a general sense of the data, and to see if it matches the graph (smallest point, largest point, average):
```
$ python3 analyzer.py program_basenames.txt
```

### ***To recreate the figure***
```
$ awk 'ORS=NR?" ":"\n"' program_basenames.txt | awk ' { print $1, $3, $6, $4, $2, $5  } ' > horizontal_basenames.txt
$ paste ngs_te_mapper_sorted_freq.txt relocate_sorted_freq.txt temp_sorted_freq.txt retroseq_sorted_freq.txt popoolatio
nte_sorted_freq.txt telocate_sorted_freq.txt > all_sorted.txt
$ cat horizontal_basenames.txt all_sorted.txt | column -t > tabular_combined.txt
$ awk ' { print $1, $2, $3, $4, $5, $6 } ' tabular_combined.txt | sed 's/ /,/g' > combined.txt
$ rm all_sorted.txt tabular_combined.txt 
```

### Time to remake the boxplot (done in RStudio)
```
> combined <- read.csv("combined.txt")
> boxplot(combined, log = 'y', las = 2, par(mar = c(7, 4.25, 0.20, 0.20)+ 0.1), col = c("#f2b57c", "#ef7d7c", "#b67eec", "#7f7cf1", "#7cb7f1", "#7af1ef"))
> mtext("Number of predicted non-reference insertions per sample", 
      side = 2, line = 3.5, cex = .9)
```
