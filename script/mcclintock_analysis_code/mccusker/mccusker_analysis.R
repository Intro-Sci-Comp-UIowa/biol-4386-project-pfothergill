# Load the required packages
require(tools)
library(xtable)
library(seqplots)

# Access the command line options:
# Argument 1 should be the path to the McCusker dataset folder
# Argument 2 should be the full path to the desired output location
args<-commandArgs(TRUE)
dir<-args[1]
output.location<-args[2]
dir.create(output.location)

ngs_te_functions<-"ngs_te_mapper_functions.R"
source(ngs_te_functions)

# Define colours for the methods
ngs_te_mapper.sr.colour<-"#EDB782"
relocate.sr.colour<-"#ED8282"
temp.sr.colour<-"#ED82B8"
temp.colour<-"#EE82EE"
temp.rp.colour<-"#B782ED"
retroseq.rp.colour<-"#8282ED"
popoolationte.rp.colour<-"#82B8ED"
telocate.rp.colour<-"#82EDED"

# List all the bed files within the of this project results
results.location<-paste(dir, "/*/*/results", sep="")
bed.files<-list.files(Sys.glob(results.location), full.names=T, pattern="nonredundant.bed$")

methods<-c("ngs_te_mapper", "relocate", "temp", "retroseq", "popoolationte", "telocate")

###########################################################################
######################## TEs at tRNA genes - table ########################
###########################################################################

# Create table of how many TEs are in or out of a window around tRNA genes

# Download and prepare tRNA annotation
system(paste("wget -q http://downloads.yeastgenome.org/curation/chromosomal_feature/archive/saccharomyces_cerevisiae_R61-1-1_20080607.gff.gz -O ", output.location, "/sacCer2annotation.gff.gz", sep=""))
system(paste("gunzip -c ", output.location, "/sacCer2annotation.gff.gz > ", output.location, "/sacCer2annotation.gff", sep=""))
system(paste("rm ", output.location, "/sacCer2annotation.gff.gz", sep=""))
system(paste("awk '{if ($3==\"tRNA\") print $0}' ", output.location, "/sacCer2annotation.gff > ", output.location, "/sacCer2tRNAs.gff", sep=""))
system(paste("rm ", output.location, "/sacCer2annotation.gff", sep=""))

methods<-c("ngs_te_mapper", "relocate", "temp", "retroseq", "popoolationte", "telocate")
TEs<-c("TY1", "TY2", "TY3", "TY3_1p", "TY4", "TY5")
tRNAtable<-vector()
system(paste("awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$4\"\\t\"$6\"\\t\"$7\"\\t\"$8\"\\t\"$9}' ", output.location, "/sacCer2tRNAs.gff > ", output.location, "/sacCer2tRNAstarts.gff", sep=""))
for(TE in TEs)
{
	TEline<-vector()
	system(paste("grep -h ", TE, "\\; ", dir, "/sacCer2/reference/reference_TE_locations_HL.gff > ", output.location, "/", TE, "carr.gff",sep=""))
	intRNA<-system(paste("bedtools window -u -l 1000 -r 500 -sw -a ", output.location, "/", TE, "carr.gff -b ", output.location, "/sacCer2tRNAstarts.gff | wc  | awk '{print $1}'", sep=""), intern=TRUE)
	total<-system(paste("wc ", output.location, "/", TE, "carr.gff | awk '{print $1}'",  sep=""), intern=TRUE)
	nottRNA<-as.numeric(total)-as.numeric(intRNA)
	#percent<-sprintf("%.2f", (as.numeric(intRNA)/as.numeric(total))*100, 2)
	percent<-round((as.numeric(intRNA)/as.numeric(total))*100)
	if (is.nan(percent)) 
	{
		percent<-"N.A."
		temethodresult<-paste(intRNA, "/", total, " (", percent, ")",sep="")
	}
	else 
	{
		temethodresult<-paste(intRNA, "/", total, " (", percent, "\\%)",sep="")
	}
	TEline<-cbind(TEline, temethodresult)
	for(method in methods)
	{
		system(paste("grep -h ", TE, "_non- ", results.location, "/*_", method, "_nonredundant.bed > ", output.location, "/", method, TE, ".list",sep=""))
		intRNA<-system(paste("bedtools window -u -l 1000 -r 500 -sw -a ", output.location, "/", method, TE, ".list -b ", output.location, "/sacCer2tRNAstarts.gff | wc  | awk '{print $1}'", sep=""), intern=TRUE)
		total<-system(paste("grep -h ", TE, "_non- ", results.location, "/*_", method, "_nonredundant.bed | wc  | awk '{print $1}'",  sep=""), intern=TRUE)
		nottRNA<-as.numeric(total)-as.numeric(intRNA)
		#percent<-sprintf("%.2f", (as.numeric(intRNA)/as.numeric(total))*100, 2)
		percent<-round((as.numeric(intRNA)/as.numeric(total))*100)
	if (is.nan(percent)) 
	{
		percent<-"N.A."
		temethodresult<-paste(intRNA, "/", total, " (", percent, ")",sep="")
	}
	else 
	{
		temethodresult<-paste(intRNA, "/", total, " (", percent, "\\%)",sep="")
	}
	TEline<-cbind(TEline, temethodresult)
		system(paste("rm ", output.location, "/", method, TE, ".list",sep=""))
	}
	tRNAtable<-rbind(tRNAtable, TEline)
}

colnames(tRNAtable)<-c("\\textbf{Carr}", "\\textbf{\\ngs}", "\\textbf{RelocaTE}", "\\textbf{TEMP}", "\\textbf{RetroSeq}", "\\textbf{PoPoolationTE}", "\\textbf{TE-locate}")
row.names(tRNAtable)<-c("\\emph{Ty1}", "\\emph{Ty2}", "\\emph{Ty3}", "\\emph{Ty3\\_1p}", "\\emph{Ty4}", "\\emph{Ty5}")

latextable1<-xtable(tRNAtable)
print(latextable1, file = paste(output.location, '/tes_at_trna_table.tex', sep=""), hline.after = c(), include.rownames = TRUE, include.colnames = FALSE, only.contents = TRUE, sanitize.text.function = function(x) x)


###########################################################################
############################ Counts box plots #############################
###########################################################################

# For each bed file
for(i in 1:length(bed.files))
{
	if (exists("dataset"))
	{	
		# If the is not the first file add the data to the previously created dataset
		temp.dataset <- read.table(bed.files[i],header=FALSE,col.names=c("Chromosome","Start","End", "TE Name", "Score", "Orientation"),skip=1,sep="\t")
		
		filename<-basename(file_path_sans_ext(bed.files[i]))
		temp.dataset$Sample<-basename(dirname(dirname(bed.files[i])))
		
		dataset<-rbind(dataset, temp.dataset)
		rm(temp.dataset)
	}
	if (!exists("dataset"))
	{
		# If this is the first file, read the data into R
		dataset <- read.table(bed.files[i],header=FALSE,col.names=c("Chromosome","Start","End", "TE Name", "Score", "Orientation"),skip=1,sep="\t")

		# Add columns for the sample name and method
		filename<-basename(file_path_sans_ext(bed.files[i]))
		dataset$Sample<-basename(dirname(dirname(bed.files[i])))
	}
}

dataset[grep("ngs_te_mapper",  dataset$TE.Name), "Method"] <- "ngs_te_mapper"
dataset[grep("relocate",  dataset$TE.Name), "Method"] <- "relocate"
dataset[grep("temp_sr",  dataset$TE.Name), "Method"] <- "tempsplitread"
dataset[grep("temp_rp",  dataset$TE.Name), "Method"] <- "tempreadpair"
dataset[grep("temp_nonab",  dataset$TE.Name), "Method"] <- "tempreference"
dataset[grep("retroseq",  dataset$TE.Name), "Method"] <- "retroseq"
dataset[grep("telocate",  dataset$TE.Name), "Method"] <- "telocate"
dataset[grep("popoolationte",  dataset$TE.Name), "Method"] <- "popoolationte"

dataset[grep("TY1",  dataset$TE.Name), "TE"] <- "Ty1"
dataset[grep("TY2",  dataset$TE.Name), "TE"] <- "Ty2"
dataset[grep("TY3",  dataset$TE.Name), "TE"] <- "Ty3"
dataset[grep("TY4",  dataset$TE.Name), "TE"] <- "Ty4"
dataset[grep("TY5",  dataset$TE.Name), "TE"] <- "Ty5"
dataset[grep("TY3_1p",  dataset$TE.Name), "TE"] <- "Ty3_1p"

# Exclude problematic sample
#no809<-subset(dataset, Sample!="SRR800809")

# Select all of the novel insertions
new.TEs<-grep("non-reference",dataset$TE.Name)

# Count the number of TEs per sample for each method
newTEs<-paste(dataset$Method[new.TEs],dataset$Sample[new.TEs])
new.TE.counts<-table(newTEs)

new.TE.counts.ngs_te_mapper<-new.TE.counts[grep("ngs_te_mapper",names(new.TE.counts))]
new.TE.counts.relocate<-new.TE.counts[grep("relocate",names(new.TE.counts))]
new.TE.counts.retroseq<-new.TE.counts[grep("retroseq",names(new.TE.counts))]
new.TE.counts.popoolationte<-new.TE.counts[grep("popoolationte",names(new.TE.counts))]
new.TE.counts.telocate<-new.TE.counts[grep("telocate",names(new.TE.counts))]
new.TE.counts.temp<-new.TE.counts[grep("tempsplitread",names(new.TE.counts))]
new.TE.counts.temprp<-new.TE.counts[grep("tempreadpair",names(new.TE.counts))]

# Plot boxplots of the counts of new TEs
pdf(paste(output.location, '/nonreflogboxplot.pdf', sep=""))
par(oma = c(0,0,0,0) + 0.1,mar = c(10,5,1,1) + 0.1)
boxplot(as.vector(new.TE.counts.ngs_te_mapper), as.vector(new.TE.counts.relocate), as.vector(new.TE.counts.temp), as.vector(new.TE.counts.temprp), as.vector(new.TE.counts.retroseq), as.vector(new.TE.counts.popoolationte), as.vector(new.TE.counts.telocate), names = c("ngs_te_mapper", "RelocaTE", "TEMP (Split-read)", "TEMP (Read-pair)", "RetroSeq", "PoPoolationTE", "TE-locate"), log = "y", las=2, at =c(1,2,3,4.5,5.5,6.5,7.5), col=c(ngs_te_mapper.sr.colour, relocate.sr.colour, temp.sr.colour, temp.rp.colour, retroseq.rp.colour, popoolationte.rp.colour, telocate.rp.colour) )
title(ylab = "Number of predicted non-reference insertions per sample", cex.lab = 1, line = 3.5)
axis(1,at=c(0.5,1,2,3,3.5),line=8.5,tick=T,labels=c("","","Split-read","",""),lwd=2,lwd.ticks=0,padj=-1.25)
axis(1,at=3+c(1,2,3,4,5),line=8.5,tick=T,labels=c("","","Read-pair","",""),lwd=2,lwd.ticks=0,padj=-1.25)
dev.off()


# Select all of the detected reference elements
old.TEs<-grep("_reference_",dataset$TE.Name)

# Count the number of TEs per sample for each method, RetroSeq is not counted as it does not detect reference TEs
oldTEs<-paste(dataset$Method[old.TEs],dataset$Sample[old.TEs])
old.TE.counts<-table(oldTEs)

old.TE.counts.ngs_te_mapper<-old.TE.counts[grep("ngs_te_mapper",names(old.TE.counts))]
old.TE.counts.relocate<-old.TE.counts[grep("relocate",names(old.TE.counts))]
old.TE.counts.temp<-old.TE.counts[grep("tempreference",names(old.TE.counts))]
old.TE.counts.popoolationte<-old.TE.counts[grep("popoolationte",names(old.TE.counts))]
old.TE.counts.telocate<-old.TE.counts[grep("telocate",names(old.TE.counts))]

# Plot boxplots of the counts of reference TEs on log scale
pdf(paste(output.location, '/reflogboxplot.pdf', sep=""))
par(oma = c(0,0,0,0) + 0.1,mar = c(9,5,1,1) + 0.1)
boxplot(as.vector(old.TE.counts.ngs_te_mapper), as.vector(old.TE.counts.relocate), as.vector(old.TE.counts.temp), as.vector(old.TE.counts.popoolationte), as.vector(old.TE.counts.telocate), names = c("ngs_te_mapper","RelocaTE", "TEMP", "PoPoolationTE", "TE-locate"), log = "y", las=2, at =c(1,2,3.5,5,6), col=c(ngs_te_mapper.sr.colour, relocate.sr.colour, temp.colour, popoolationte.rp.colour, telocate.rp.colour))
title(ylab = "Number of predicted reference TEs per sample", cex.lab = 1, line = 3.5)
axis(1,at=c(0.5,1.5,2.5),line=7.5,tick=T,labels=c("","Split-read",""),lwd=2,lwd.ticks=0,padj=-1.25)
axis(1,at=c(2.8,3.5,4.2),line=7.5,tick=T,labels=c("","Non-absent",""),lwd=2,lwd.ticks=0,padj=-1.25)
axis(1,at=4.5+c(0,1,2),line=7.5,tick=T,labels=c("","Read-pair",""),lwd=2,lwd.ticks=0,padj=-1.25)
dev.off()

###########################################################################
############################# TSD lengths #################################
###########################################################################

# Create sequence logos of TSDs and insertion sites
# Use GetFasta from ngs_te_mapper to import reference sequence
reference.location<-paste(dir, "/sacCer2/reference/sacCer2.fasta", sep="")
reference<-GetFasta(reference.location)

# Select only the new insertions for ngs_te_mapper and RelocaTE (base pair accurate insertions)
ngs_te_mapper.insertions<-dataset[dataset$Method == "ngs_te_mapper",]
new.ngs_te_mapper.insertions<-ngs_te_mapper.insertions[grep("non-reference", ngs_te_mapper.insertions$TE.Name),]

relocate.insertions<-dataset[dataset$Method == "relocate",]
new.relocate.insertions<-relocate.insertions[grep("non-reference", relocate.insertions$TE.Name),]

# Select SR insertions from TEMP - have TSD - Or just select with correct TSD length later?
temp.insertions<-dataset[dataset$Method == "tempsplitread",]
new.temp.insertions<-temp.insertions[grep("non-reference", temp.insertions$TE.Name),]

# Add a TSD calculation to all data
new.temp.insertions$tsd<-new.temp.insertions$End-temp.insertions$Start
new.relocate.insertions$tsd<-new.relocate.insertions$End-new.relocate.insertions$Start
new.ngs_te_mapper.insertions$tsd<-new.ngs_te_mapper.insertions$End-new.ngs_te_mapper.insertions$Start

# Get only unique instances of TE insertions
uniq.temp.insertions<-new.temp.insertions[!duplicated(temp.insertions[c("Chromosome","Start","End","TE")]),]
uniq.ngs_te_mapper.insertions<-new.ngs_te_mapper.insertions[!duplicated(new.ngs_te_mapper.insertions[c("Chromosome","Start","End","TE")]),]
uniq.relocate.insertions<-new.relocate.insertions[!duplicated(new.relocate.insertions[c("Chromosome","Start","End","TE")]),]

# Select only new TEs, not including TY3_1p (can not transpose) or TY5 (no results or these methods) to analyse
TE.Names<-names(table(dataset$TE))
TE.Names<-TE.Names[grep("1p", TE.Names, invert = TRUE)]
TE.Names<-TE.Names[grep("Ty5", TE.Names, invert = TRUE)]

pdf(paste(output.location, '/TSDlengths.pdf', sep=""))
par(mfcol=c(4,3), mar=c(2.5,2,1,1), oma=c(4,2.5,2,0))

# ngs_te_mapper TSD lengths
for ( i in 1:length(TE.Names))
{
	TE.tsds<-which(uniq.ngs_te_mapper.insertions$TE == TE.Names[i])
	uniq.ngs_te_mapper.insertions.matrix<-matrix(uniq.ngs_te_mapper.insertions[TE.tsds,10])
	if(TE.Names[i] == "Ty1")
	{
		hist(uniq.ngs_te_mapper.insertions.matrix, xlab="", ylab="Number of insertions", xlim=c(0,15), ylim=c(0,40), main="ngs_te_mapper", breaks=seq(min(uniq.ngs_te_mapper.insertions.matrix)-0.5, max(uniq.ngs_te_mapper.insertions.matrix)+0.5, by=1), col=c(ngs_te_mapper.sr.colour))
	}
	if(TE.Names[i] == "Ty2")
	{
		hist(uniq.ngs_te_mapper.insertions.matrix, xlab="", ylab="Number of insertions", xlim=c(0,15), ylim=c(0,170), main="", breaks=seq(min(uniq.ngs_te_mapper.insertions.matrix)-0.5, max(uniq.ngs_te_mapper.insertions.matrix)+0.5, by=1), col=c(ngs_te_mapper.sr.colour))
	}
	if(TE.Names[i] == "Ty3")
	{
			hist(uniq.ngs_te_mapper.insertions.matrix, xlab="", ylab="Number of insertions", xlim=c(0,15), ylim=c(0,75), main="", breaks=seq(min(uniq.ngs_te_mapper.insertions.matrix)-0.5, max(uniq.ngs_te_mapper.insertions.matrix)+0.5, by=1), col=c(ngs_te_mapper.sr.colour))
	}
	if(TE.Names[i] == "Ty4")
	{
			hist(uniq.ngs_te_mapper.insertions.matrix, xlab="TSD length (bp)", ylab="Number of insertions", xlim=c(0,15), main="", ylim=c(0,80), breaks=seq(min(uniq.ngs_te_mapper.insertions.matrix)-0.5, max(uniq.ngs_te_mapper.insertions.matrix)+0.5, by=1), col=c(ngs_te_mapper.sr.colour))
	}
}

# RelocaTE TSD lengths
for ( i in 1:length(TE.Names))
{
	TE.tsds<-which(uniq.relocate.insertions$TE == TE.Names[i])
	uniq.relocate.insertions.matrix<-matrix(uniq.relocate.insertions[TE.tsds,10])
	if(TE.Names[i] == "Ty1")
	{
		hist(uniq.relocate.insertions.matrix, xlab="",ylab="", xlim=c(0,15), ylim=c(0,40), main="RelocaTE", breaks=seq(min(uniq.relocate.insertions.matrix)-0.5, max(uniq.relocate.insertions.matrix)+0.5, by=1), col=c(relocate.sr.colour))
	}
	if(TE.Names[i] == "Ty2")
	{
		hist(uniq.relocate.insertions.matrix, xlab="",ylab="", xlim=c(0,15), ylim=c(0,170), main="", breaks=seq(min(uniq.relocate.insertions.matrix)-0.5, max(uniq.relocate.insertions.matrix)+0.5, by=1), col=c(relocate.sr.colour))
	}
	if(TE.Names[i] == "Ty3")
	{
		hist(uniq.relocate.insertions.matrix, xlab="",ylab="", xlim=c(0,15), ylim=c(0,75), main="", breaks=seq(min(uniq.relocate.insertions.matrix)-0.5, max(uniq.relocate.insertions.matrix)+0.5, by=1), col=c(relocate.sr.colour))
	}
	if(TE.Names[i] == "Ty4")
	{
		hist(uniq.relocate.insertions.matrix, xlab="TSD length (bp)",ylab="", xlim=c(0,15), ylim=c(0,80), main="", breaks=seq(min(uniq.relocate.insertions.matrix)-0.5, max(uniq.relocate.insertions.matrix)+0.5, by=1), col=c(relocate.sr.colour))
	}
}

# TEMP TSD lengths
for ( i in 1:length(TE.Names))
{
	TE.tsds<-which(uniq.temp.insertions$TE == TE.Names[i])
	uniq.temp.insertions.matrix<-matrix(uniq.temp.insertions[TE.tsds,10])
	if(TE.Names[i] == "Ty1")
	{
		hist(uniq.temp.insertions.matrix, xlab="", ylab="", xlim=c(0,15), ylim=c(0,40),  main="TEMP (SR)", breaks=seq(min(uniq.temp.insertions.matrix)-0.5, max(uniq.temp.insertions.matrix)+0.5, by=1), col=c(temp.sr.colour))
	}
	if(TE.Names[i] == "Ty2")
	{
		hist(0, xlab="", ylab="", xlim=c(0,15), ylim=c(0,170), col=c(temp.sr.colour), main="")
	}
	if(TE.Names[i] == "Ty3")
	{
		hist(uniq.temp.insertions.matrix, xlab="", ylab="", xlim=c(0,15), ylim=c(0,75), main="", breaks=seq(min(uniq.temp.insertions.matrix)-0.5, max(uniq.temp.insertions.matrix)+0.5, by=1), col=c(temp.sr.colour))
	}
	if(TE.Names[i] == "Ty4")
	{
		hist(uniq.temp.insertions.matrix, xlab="TSD length (bp)", ylab="", xlim=c(0,15), ylim=c(0,80), main="", breaks=seq(min(uniq.temp.insertions.matrix)-0.5, max(uniq.temp.insertions.matrix)+0.5, by=1), col=c(temp.sr.colour))
	}
}

mtext(substitute(paste(italic('Ty4'))), side = 2, line = 0.5, outer = TRUE, at=0.14)
mtext(substitute(paste(italic('Ty3'))), side = 2, line = 0.5, outer = TRUE, at=0.37)
mtext(substitute(paste(italic('Ty2'))), side = 2, line = 0.5, outer = TRUE, at=0.63)
mtext(substitute(paste(italic('Ty1'))), side = 2, line = 0.5, outer = TRUE, at=0.89)

mtext("Predicted TSD lengths (bp)", side = 1, line = 0.5, outer = TRUE)

dev.off()


###########################################################################
##################### tRNA correlation seqplots ###########################
###########################################################################

# Do Seqplots of insertions around tRNA genes

# Create files for seqplots - one bigwig per method per TE family (with results?)
srmethods<-c("ngs_te_mapper", "relocate", "temp")
TEs<-c("TY1", "TY2", "TY3", "TY4")
system(paste("cut -f1,2 ", dir, "/sacCer2/reference/sacCer2.fasta.fai > ", output.location, "/sacCer2.lengths", sep=""))

pdf(paste(output.location, '/tRNA_sr_correlation.pdf', sep=""), width=8)
par(mfcol=c(4,3), mar=c(1.2,1.2,1,1), oma=c(2.1,2.5,0.5,0))

for (method in srmethods)
{
	for (TE in TEs)
	{
		if ((TE  == "TY2") && (method == "temp"))
		{
			system(paste("awk '{print $1\"\t0\t1\t0\"}' ", output.location, "/sacCer2.lengths > ", output.location, "/", TE, "_sr_", method, ".bed", sep=""))
		}
		else
		{
			system(paste("grep --no-filename ", TE, "_non- ", results.location, "/*", method, "_nonredundant.bed | grep _sr_ | sort -k 1,1 > ", output.location, "/", TE, "_sr_", method, ".bed", sep=""))
		}
		system(paste("bedtools genomecov -i ", output.location, "/", TE, "_sr_", method, ".bed -g ", output.location, "/sacCer2.lengths -bga -trackline | wigToBigWig -clip stdin ", output.location, "/sacCer2.lengths ", output.location, "/", TE, "_sr_", method, ".bw", sep=""))
		#system(paste("rm ", output.location, "/", TE, "_sr_", method, ".bed", sep=""))

		plotSet<-getPlotSetArray(paste(output.location, "/", TE, "_sr_", method, ".bw", sep=""), paste(output.location, "/sacCer2tRNAs.gff", sep=""), 'sacCer2', bin = 1L, xmin = 1000L, xmax = 500L, type="pf", add_heatmap=FALSE)

		if ((TE  == "TY2") && (method == "temp"))
		{
			plotAverage(plotSet, main="", legend=FALSE, ylim=c(0,1), cex.main=10, cex.axis=10, error.estimates=FALSE, colvec=get(paste(method, ".sr.colour", sep="")))
		}
		else
		{
			plotAverage(plotSet, main="", legend=FALSE, ylim=c(0,max(plotSet$data[[1]][[1]]$means)*1.05), cex.main=10, cex.axis=10, error.estimates=FALSE, colvec=get(paste(method, ".sr.colour", sep="")))
		}
	}
}

mtext("Coordinates relative to tRNA start", side = 1, line = 0.7, outer = TRUE)
mtext(substitute(paste(italic('Ty4'))), side = 2, line = 0.5, outer = TRUE, at=0.15)
mtext(substitute(paste(italic('Ty3'))), side = 2, line = 0.5, outer = TRUE, at=0.40)
mtext(substitute(paste(italic('Ty2'))), side = 2, line = 0.5, outer = TRUE, at=0.63)
mtext(substitute(paste(italic('Ty1'))), side = 2, line = 0.5, outer = TRUE, at=0.88)
mtext("ngs_te_mapper", side = 3, line = -0.75, outer = TRUE, at=0.20)
mtext("RelocaTE", side = 3, line = -0.75, outer = TRUE, at=0.53)
mtext("TEMP (SR)", side = 3, line = -0.75, outer = TRUE, at=0.85)

dev.off()

rpmethods<-c("temp", "retroseq", "popoolationte", "telocate")
pdf(paste(output.location, '/tRNA_rp_correlation.pdf', sep=""), width=11.34, height=7)
par(mfcol=c(4,4), mar=c(1.2,1.2,1,1), oma=c(2.1,2.5,0.5,0))

for (method in rpmethods)
{
	for (TE in TEs)
	{
		system(paste("grep --no-filename ", TE, "_non- ", results.location, "/*", method, "_nonredundant.bed | grep _rp_ | sort -k 1,1 > ", output.location, "/", TE, "_rp_", method, ".bed", sep=""))
		system(paste("bedtools genomecov -i ", output.location, "/", TE, "_rp_", method, ".bed -g ", output.location, "/sacCer2.lengths -bga -trackline | wigToBigWig -clip stdin ", output.location, "/sacCer2.lengths ", output.location, "/", TE, "_rp_", method, ".bw", sep=""))
		system(paste("rm ", output.location, "/", TE, "_rp_", method, ".bed", sep=""))
		plotSet<-getPlotSetArray(paste(output.location, "/", TE, "_rp_", method, ".bw", sep=""), paste(output.location, "/sacCer2tRNAs.gff", sep=""), 'sacCer2', bin = 1L, xmin = 1000L, xmax = 500L, type="pf", add_heatmap=FALSE)
		plotAverage(plotSet, main="", legend=FALSE, ylim=c(0,max(plotSet$data[[1]][[1]]$means)*1.05), cex.main=10, cex.axis=10, error.estimates=FALSE, colvec=get(paste(method, ".rp.colour", sep="")))
	}
}

mtext("Coordinates relative to tRNA start", side = 1, line = 0.7, outer = TRUE)
mtext(substitute(paste(italic('Ty4'))), side = 2, line = 0.5, outer = TRUE, at=0.15)
mtext(substitute(paste(italic('Ty3'))), side = 2, line = 0.5, outer = TRUE, at=0.40)
mtext(substitute(paste(italic('Ty2'))), side = 2, line = 0.5, outer = TRUE, at=0.63)
mtext(substitute(paste(italic('Ty1'))), side = 2, line = 0.5, outer = TRUE, at=0.88)
mtext("TEMP (RP)", side = 3, line = -0.75, outer = TRUE, at=0.15)
mtext("RetroSeq", side = 3, line = -0.75, outer = TRUE, at=0.39)
mtext("PoPoolationTE", side = 3, line = -0.75, outer = TRUE, at=0.64)
mtext("TE-locate", side = 3, line = -0.75, outer = TRUE, at=0.89)

dev.off()

################## Save workspace and exit the script #####################
save.image(file=paste(output.location, "/workspace.Rdata", sep=""))
