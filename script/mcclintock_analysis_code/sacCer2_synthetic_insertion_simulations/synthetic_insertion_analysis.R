# Load the required packages
require(tools)
library(xtable)
library(seqplots)
library(VennDiagram)
library(grid)
library(gridBase)
library(lattice)
library(gridExtra)
options(error=traceback)

# Argument 1 should be the path to the singleinsertionsims folder
# Argument 2 should be the full path to the desired output location
args<-commandArgs(TRUE)
dir<-args[1]
output.location<-args[2]

ngs_te_functions<-paste("ngs_te_mapper_functions.R")
source(ngs_te_functions)

# Define colours for the methods
ngs_te_mapper.sr.colour<-"#EDB78299"
ngs_te_mapper.sr.colourpos<-"#D49E6999"
ngs_te_mapper.sr.colourneg<-"#FFD19C99"
relocate.sr.colour<-"#ED828299"
relocate.sr.colourpos<-"#D4696999"
relocate.sr.colourneg<-"#FF9C9C99"
temp.sr.colour<-"#ED82B899"
temp.sr.colourpos<-"#D4699F99"
temp.sr.colourneg<-"#FF9CD299"
temp.colour<-"#EE82EE99"
temp.rp.colour<-"#B782ED99"
temp.rp.colourpos<-"#9E69D499"
temp.rp.colourneg<-"#D19CFF8C"
retroseq.rp.colour<-"#8282ED99"
retroseq.rp.colourpos<-"#6969D499"
retroseq.rp.colourneg<-"#9C9CFF8C"
popoolationte.rp.colour<-"#82B8ED99"
popoolationte.rp.colourpos<-"#699FD499"
popoolationte.rp.colourneg<-"#9CD2FF8C"
telocate.rp.colour<-"#82EDED99"
telocate.rp.colourpos<-"#69D4D499"
telocate.rp.colourneg<-"#9CFFFF8C"

# Create an output directory for results
dir.create(output.location)
results.location<-paste(dir, "/*/*/results", sep="")

reverse.dir<-paste(dir, "/../singlerevinsertionsims", sep="")
reverse.results.location<-paste(reverse.dir, "/*/*/results", sep="")

###########################################################################
################### Split read single insertion accuracy ##################
###########################################################################

# Create files for seqplots - one bigwig per method per TE family (with results?)
srmethods<-c("ngs_te_mapper", "relocate", "temp")
TEs<-c("TY1", "TY2", "TY3", "TY4")
system(paste("cut -f1,2 ", dir, "/sacCer2/reference/sacCer2.fasta.fai > ", output.location, "/sacCer2.lengths", sep=""))

# Prepare records of tRNA locations (add one to the start because genomecov in 1based coords)
system(paste("awk '{if ($2==0) printf $1\"\t\"$3\"\t\"; else print $2+1}' ", dir, "/data/insertions.bed | awk 'NR%4==1{print $1\"\t\"$3\"\t\"$2\"\tTy1\"}NR%4==2{print $1\"\t\"$3\"\t\"$2\"\tTy2\"}NR%4==3{print $1\"\t\"$3\"\t\"$2\"\tTy3\"}NR%4==0{print $1\"\t\"$3\"\t\"$2\"\tTy4\"}' >", output.location, "/singleinsertionlocations.bed", sep=""))

system(paste("grep Ty1 ", output.location, "/singleinsertionlocations.bed > ", output.location, "/singleinsertionlocationsTY1.bed", sep=""))
system(paste("grep Ty2 ", output.location, "/singleinsertionlocations.bed > ", output.location, "/singleinsertionlocationsTY2.bed", sep=""))
system(paste("grep Ty3 ", output.location, "/singleinsertionlocations.bed > ", output.location, "/singleinsertionlocationsTY3.bed", sep=""))
system(paste("grep Ty4 ", output.location, "/singleinsertionlocations.bed > ", output.location, "/singleinsertionlocationsTY4.bed", sep=""))

pdf(paste(output.location, '/singleins_sr_correlation.pdf', sep=""))
par(mfcol=c(4,3), mar=c(1.2,1.2,1,1), oma=c(2.4,2.4,0.5,0))

for (method in srmethods)
{
	for (TE in TEs)
	{
		if ((TE  == "TY2") && (method == "temp"))
		{
			system(paste("awk '{print $1\"\t0\t1\t0\"}' ", output.location, "/sacCer2.lengths > ", output.location, "/", TE, "+_sr_", method, ".bed", sep=""))
			system(paste("awk '{print $1\"\t0\t1\t0\"}' ", output.location, "/sacCer2.lengths > ", output.location, "/", TE, "-_sr_", method, ".bed", sep=""))
		}
		else
		{
			system(paste("grep --no-filename ", TE, "_non- ", results.location, "/*", method, "_nonredundant.bed | grep _sr_ | sort -k 1,1 > ", output.location, "/", TE, "+_sr_", method, ".bed", sep=""))
			system(paste("grep --no-filename ", TE, "_non- ", reverse.results.location, "/*", method, "_nonredundant.bed | grep _sr_ | sort -k 1,1 > ", output.location, "/", TE, "-_sr_", method, ".bed", sep=""))
		}
		system(paste("bedtools genomecov -i ", output.location, "/", TE, "+_sr_", method, ".bed -g ", output.location, "/sacCer2.lengths -bga -trackline >", output.location, "/", TE, "+_sr_", method, ".bedgraph", sep=""))
		system(paste("bedtools genomecov -i ", output.location, "/", TE, "-_sr_", method, ".bed -g ", output.location, "/sacCer2.lengths -bga -trackline >", output.location, "/", TE, "-_sr_", method, ".bedgraph", sep="")) 

		system(paste("awk '{if($4>1) print $1\"\t\"$2\"\t\"$3\"\t0\"; else print $0}' ", output.location, "/", TE, "+_sr_", method, ".bedgraph >", output.location, "/", TE, "+_sr_", method, "edited.bedgraph", sep=""))	
		system(paste("awk '{if($4>1) print $1\"\t\"$2\"\t\"$3\"\t0\"; else print $0}' ", output.location, "/", TE, "-_sr_", method, ".bedgraph >", output.location, "/", TE, "-_sr_", method, "edited.bedgraph", sep=""))
		
		print(system(paste("awk '{if($4>1) print $1\"\t\"$2\"\t\"$3\"\t\"$4}' ", output.location, "/", TE, "+_sr_", method, ".bedgraph", sep="")))
		print(system(paste("awk '{if($4>1) print $1\"\t\"$2\"\t\"$3\"\t\"$4}' ", output.location, "/", TE, "-_sr_", method, ".bedgraph", sep="")))
			
		system(paste("wigToBigWig -clip ", output.location, "/", TE, "+_sr_", method, "edited.bedgraph ", output.location, "/sacCer2.lengths ", output.location, "/", TE, "+_sr_", method, ".bw", sep=""))
		system(paste("wigToBigWig -clip ", output.location, "/", TE, "-_sr_", method, "edited.bedgraph ", output.location, "/sacCer2.lengths ", output.location, "/", TE, "-_sr_", method, ".bw", sep=""))

		system(paste("rm ", output.location, "/", TE, "+_sr_", method, ".bed ", output.location, "/", TE, "-_sr_", method, ".bed", sep=""))
		system(paste("rm ", output.location, "/", TE, "+_sr_", method, ".bedgraph ", output.location, "/", TE, "-_sr_", method, ".bedgraph", sep=""))
		system(paste("rm ", output.location, "/", TE, "+_sr_", method, "edited.bedgraph ", output.location, "/", TE, "-_sr_", method, "edited.bedgraph", sep=""))
	
		bwlist<-c(paste(output.location, "/", TE, "-_sr_", method, ".bw", sep=""),paste(output.location, "/", TE, "+_sr_", method, ".bw", sep="")) 
		plotSetpos<-getPlotSetArray(bwlist, paste(output.location, "/singleinsertionlocations", TE, ".bed", sep=""), 'sacCer2', bin = 1L, xmin = 10L, xmax = 15L, type="pf", add_heatmap=FALSE)
		listy<-unlist(plotSetpos)
		if ((TE  == "TY2") && (method == "temp"))
		{
			plotAverage(listy, legend=FALSE, ylim=c(0,1), cex.main=10, cex.axis=10,  error.estimates=FALSE, colvec=get(paste(method, ".sr.colour", sep="")))
		}
		else
		{
			plotAverage(listy, legend=FALSE, ylim=c(0,1), cex.main=10, cex.axis=10,  error.estimates=FALSE, colvec=c(get(paste(method, ".sr.colourneg", sep="")), get(paste(method, ".sr.colourpos", sep=""))))
		}
	}
}

mtext("Coordinates relative to start of TSD produced by a simulated insertion", side = 1, line = 0.7, outer = TRUE)
mtext(substitute(paste(italic('Ty4'))), side = 2, line = 0.5, outer = TRUE, at=0.14)
mtext(substitute(paste(italic('Ty3'))), side = 2, line = 0.5, outer = TRUE, at=0.38)
mtext(substitute(paste(italic('Ty2'))), side = 2, line = 0.5, outer = TRUE, at=0.63)
mtext(substitute(paste(italic('Ty1'))), side = 2, line = 0.5, outer = TRUE, at=0.87)
mtext("ngs_te_mapper", side = 3, line = -0.75, outer = TRUE, at=0.20)
mtext("RelocaTE", side = 3, line = -0.75, outer = TRUE, at=0.53)
mtext("TEMP (SR)", side = 3, line = -0.75, outer = TRUE, at=0.85)

dev.off()

#############################################################################
##################### Read pair single insertion accuracy ###################
#############################################################################

rpmethods<-c("temp", "retroseq", "popoolationte", "telocate")
pdf(paste(output.location, '/singleins_rp_correlation.pdf', sep=""), width=11.34, height=7)
par(mfcol=c(4,4), mar=c(1.2,1.2,1,1), oma=c(2.4,2.4,0.5,0))

for (method in rpmethods)
{
	count<-1
	plotSets<-list()
	for (TE in TEs)
	{
		system(paste("grep --no-filename ", TE, "_non- ", results.location, "/*", method, "_nonredundant.bed | grep _rp_ | sort -k 1,1 > ", output.location, "/", TE, "+_rp_", method, ".bed", sep=""))
		system(paste("grep --no-filename ", TE, "_non- ", reverse.results.location, "/*", method, "_nonredundant.bed | grep _rp_ | sort -k 1,1 > ", output.location, "/", TE, "-_rp_", method, ".bed", sep=""))

		system(paste("bedtools genomecov -i ", output.location, "/", TE, "+_rp_", method, ".bed -g ", output.location, "/sacCer2.lengths -bga -trackline >", output.location, "/", TE, "+_rp_", method, ".bedgraph", sep=""))
		system(paste("bedtools genomecov -i ", output.location, "/", TE, "-_rp_", method, ".bed -g ", output.location, "/sacCer2.lengths -bga -trackline >", output.location, "/", TE, "-_rp_", method, ".bedgraph", sep=""))

		system(paste("awk '{if($4>1) print $1\"\t\"$2\"\t\"$3\"\t0\"; else print $0}' ", output.location, "/", TE, "+_rp_", method, ".bedgraph >", output.location, "/", TE, "+_rp_", method, "edited.bedgraph", sep=""))
		system(paste("awk '{if($4>1) print $1\"\t\"$2\"\t\"$3\"\t0\"; else print $0}' ", output.location, "/", TE, "-_rp_", method, ".bedgraph >", output.location, "/", TE, "-_rp_", method, "edited.bedgraph", sep=""))

		print(system(paste("awk '{if($4>1) print $1\"\t\"$2\"\t\"$3\"\t\"$4}' ", output.location, "/", TE, "+_rp_", method, ".bedgraph", sep="")))
		print(system(paste("awk '{if($4>1) print $1\"\t\"$2\"\t\"$3\"\t\"$4}' ", output.location, "/", TE, "-_rp_", method, ".bedgraph", sep="")))

		system(paste("wigToBigWig -clip ", output.location, "/", TE, "+_rp_", method, "edited.bedgraph ", output.location, "/sacCer2.lengths ", output.location, "/", TE, "+_rp_", method, ".bw", sep=""))
		system(paste("wigToBigWig -clip ", output.location, "/", TE, "-_rp_", method, "edited.bedgraph ", output.location, "/sacCer2.lengths ", output.location, "/", TE, "-_rp_", method, ".bw", sep=""))

		bwlist<-c(paste(output.location, "/", TE, "+_rp_", method, ".bw", sep=""), paste(output.location, "/", TE, "-_rp_", method, ".bw", sep=""))
		plotSet<-getPlotSetArray(bwlist, paste(output.location, "/singleinsertionlocations", TE, ".bed", sep=""), 'sacCer2', bin = 1L, xmin = 500L, xmax = 505L, type="pf", add_heatmap=FALSE)
		plotSets[count]<-unlist(plotSet)
		count<-count+1
	}
	maxscale<-0
	for(count in 1:4)
	{
		tempmax<-max(plotSets[[1]]$data[[1]][[1]]$means, plotSets[[1]]$data[[1]][[2]]$means)
		maxscale<-max(maxscale, tempmax)
	}
	maxscale<-maxscale*1.05
	count<-1
	for (TE in TEs)
	{
		if (method == "popoolationte")
		{
			plotAverage(plotSets[[count]], main="", legend=FALSE, ylim=c(0, 1), cex.main=10, cex.axis=10, error.estimates=FALSE, colvec=c(get(paste(method, ".rp.colourpos", sep="")), get(paste(method, ".rp.colourneg", sep=""))))
		}
		else
		{

			plotAverage(plotSets[[count]], main="", legend=FALSE, ylim=c(0, maxscale), cex.main=10, cex.axis=10, error.estimates=FALSE, colvec=c(get(paste(method, ".rp.colourpos", sep="")), get(paste(method, ".rp.colourneg", sep=""))))
		}	
		count<-count+1
	}
}

mtext("Coordinates relative to start of TSD produced by a simulated insertion", side = 1, line = 0.7, outer = TRUE)
mtext(substitute(paste(italic('Ty4'))), side = 2, line = 0.5, outer = TRUE, at=0.14)
mtext(substitute(paste(italic('Ty3'))), side = 2, line = 0.5, outer = TRUE, at=0.39)
mtext(substitute(paste(italic('Ty2'))), side = 2, line = 0.5, outer = TRUE, at=0.63)
mtext(substitute(paste(italic('Ty1'))), side = 2, line = 0.5, outer = TRUE, at=0.89)
mtext("TEMP (RP)", side = 3, line = -0.75, outer = TRUE, at=0.15)
mtext("RetroSeq", side = 3, line = -0.75, outer = TRUE, at=0.39)
mtext("PoPoolationTE", side = 3, line = -0.75, outer = TRUE, at=0.64)
mtext("TE-locate", side = 3, line = -0.75, outer = TRUE, at=0.89)

dev.off()

###########################################################################
################### Single Insertion Simulations ##########################
###########################################################################

nonreftable<-matrix(, nrow = 6, ncol = 13)

nonreftable[1,1]<-"Reference TEs mean"
nonreftable[2,1]<-"Non-reference TEs mean"

methods=array(c("ngs_te_mapper", "relocate", "temp", "retroseq", "popoolationte", "telocate"))
thisnonrefrowmean<-vector()
for (method in methods)
{
	nonrefmethodcount<-vector()
	nonrefmethodcountrev<-vector()
	for(rep in 1:299)
	{
		nonrefmethodcount<-cbind(nonrefmethodcount, system(paste("grep non-reference ", dir, "/sacCer2/insertion", rep, "/results/insertion", rep, "_", method, "_nonredundant.bed | wc | awk '{print $1}'", sep=""), intern=TRUE))
		nonrefmethodcountrev<-cbind(nonrefmethodcountrev, system(paste("grep non-reference ", reverse.dir, "/sacCer2/insertion", rep, "/results/insertion", rep, "_", method, "_nonredundant.bed | wc | awk '{print $1}'", sep=""), intern=TRUE))
	}
	nonrefmethodcountresmean<-c(sprintf("%.2f", round(mean(as.numeric(nonrefmethodcount)), 2)), sprintf("%.2f", round(mean(as.numeric(nonrefmethodcountrev)), 2)))
	thisnonrefrowmean<-cbind(thisnonrefrowmean, nonrefmethodcountresmean)
}

thisrefrowmean<-vector()
for (method in methods)
{
	refmethodcount<-vector()
	refmethodcountrev<-vector()
	for(rep in 1:299)
	{
		refmethodcount<-cbind(refmethodcount, system(paste("grep _reference ", dir, "/sacCer2/insertion", rep, "/results/insertion", rep, "_", method, "_nonredundant.bed | wc | awk '{print $1}'", sep=""), intern=TRUE))
		refmethodcountrev<-cbind(refmethodcount, system(paste("grep _reference ", reverse.dir, "/sacCer2/insertion", rep, "/results/insertion", rep, "_", method, "_nonredundant.bed | wc | awk '{print $1}'", sep=""), intern=TRUE))
	}
	thisrefrowmean<-cbind(thisrefrowmean, sprintf("%.2f", round(mean(as.numeric(refmethodcount)), 2)))
	thisrefrowmean<-cbind(thisrefrowmean, sprintf("%.2f", round(mean(as.numeric(refmethodcountrev)), 2)))
}

nonreftable[1,2:13]<-thisrefrowmean
nonreftable[2,2:13]<-thisnonrefrowmean
nonreftable[1,8]<-"N.A."
nonreftable[1,9]<-"N.A."

nonreftable[3,1]<-"Exact"
nonreftable[4,1]<-"Within 100bp"
nonreftable[5,1]<-"Within 300bp"
nonreftable[6,1]<-"Within 500bp"

# Create unedited insertions locations
system(paste("awk '{if ($2==0) printf $1\"\t\"$3\"\t\"; else print $2}' ", dir, "/data/insertions.bed | awk 'NR%4==1{print $1\"\t\"$3\"\t\"$2\"\tTY1\"}NR%4==2{print $1\"\t\"$3\"\t\"$2\"\tTY2\"}NR%4==3{print $1\"\t\"$3\"\t\"$2\"\tTY3\"}NR%4==0{print $1\"\t\"$3\"\t\"$2\"\tTY4\"}' > ", output.location, "/singleinsertionlocations.bed", sep=""))

m<-2
for (method in methods)
{
	exactcount<-0
	within100count<-0
	within300count<-0
	within500count<-0
	exactcountneg<-0
	within100countneg<-0
	within300countneg<-0
	within500countneg<-0
	for (rep in 1:299)
	{
		system(paste("sed '", rep, "q;d' " , output.location, "/singleinsertionlocations.bed > ", output.location, "/line.bed", sep="")
)
		exactcount<-exactcount + as.numeric(system(paste("awk -F '_non' '{print $1}' ", dir, "/sacCer2/insertion", rep, "/results/insertion", rep, "_", method, "_nonredundant.bed | sort | comm -12 - ", output.location, "/line.bed | wc | awk '{print $1}'", sep=""), intern=TRUE))
		exactcountneg<-exactcountneg + as.numeric(system(paste("awk -F '_non' '{print $1}' ", reverse.dir, "/sacCer2/insertion", rep, "/results/insertion", rep, "_", method, "_nonredundant.bed | sort | comm -12 - ", output.location, "/line.bed | wc | awk '{print $1}'", sep=""), intern=TRUE))

		within100count<-within100count + as.numeric(system(paste("bedtools window -w 100 -a ", output.location, "/line.bed -b ", dir, "/sacCer2/insertion", rep, "/results/insertion", rep, "_", method, "_nonredundant.bed | awk -F'_non' '{print $1}' | awk '{if ($4 == $8) print $0}' | wc | awk '{print $1}'", sep=""), intern=TRUE))
		within100countneg<-within100countneg + as.numeric(system(paste("bedtools window -w 100 -a ", output.location, "/line.bed -b ", reverse.dir, "/sacCer2/insertion", rep, "/results/insertion", rep, "_", method, "_nonredundant.bed | awk -F'_non' '{print $1}' | awk '{if ($4 == $8) print $0}' | wc | awk '{print $1}'", sep=""), intern=TRUE))

		within300count<-within300count + as.numeric(system(paste("bedtools window -w 300 -a ", output.location, "/line.bed -b ", dir, "/sacCer2/insertion", rep, "/results/insertion", rep, "_", method, "_nonredundant.bed | awk -F'_non' '{print $1}' | awk '{if ($4 == $8) print $0}' | wc | awk '{print $1}'", sep=""), intern=TRUE))
		within300countneg<-within300countneg + as.numeric(system(paste("bedtools window -w 300 -a ", output.location, "/line.bed -b ", reverse.dir, "/sacCer2/insertion", rep, "/results/insertion", rep, "_", method, "_nonredundant.bed | awk -F'_non' '{print $1}' | awk '{if ($4 == $8) print $0}' | wc | awk '{print $1}'", sep=""), intern=TRUE))

		within500count<-within500count + as.numeric(system(paste("bedtools window -w 500 -a ", output.location, "/line.bed -b ", dir, "/sacCer2/insertion", rep, "/results/insertion", rep, "_", method, "_nonredundant.bed | awk -F'_non' '{print $1}' | awk '{if ($4 == $8) print $0}' | wc | awk '{print $1}'", sep=""), intern=TRUE))
		within500countneg<-within500countneg + as.numeric(system(paste("bedtools window -w 500 -a ", output.location, "/line.bed -b ", reverse.dir, "/sacCer2/insertion", rep, "/results/insertion", rep, "_", method, "_nonredundant.bed | awk -F'_non' '{print $1}' | awk '{if ($4 == $8) print $0}' | wc | awk '{print $1}'", sep=""), intern=TRUE))

	}
	print(paste(method, " ", exactcount))
	propexact<-exactcount / 299
	prop100<-within100count / 299
	prop300<-within300count / 299
	prop500<-within500count / 299
	propexactneg<-exactcountneg / 299
	prop100neg<-within100countneg / 299
	prop300neg<-within300countneg / 299
	prop500neg<-within500countneg / 299
	nonreftable[3,m]<-sprintf("%.2f", round(propexact, 2))
	nonreftable[3,m+1]<-sprintf("%.2f", round(propexactneg, 2))
	nonreftable[4,m]<-sprintf("%.2f", round(prop100, 2))
	nonreftable[4,m+1]<-sprintf("%.2f", round(prop100neg, 2))
	nonreftable[5,m]<-sprintf("%.2f", round(prop300, 2))
	nonreftable[5,m+1]<-sprintf("%.2f", round(prop300neg, 2))
	nonreftable[6,m]<-sprintf("%.2f", round(prop500, 2))
	nonreftable[6,m+1]<-sprintf("%.2f", round(prop500neg, 2))
	m<-m+2
}

latextable1<-xtable(nonreftable)
print(latextable1, file = paste(output.location, '/singleins_table.tex', sep=""), hline.after = c(), include.rownames = FALSE, include.colnames = FALSE, only.contents = TRUE, sanitize.text.function = function(x) x)

###########################################################################
################################ TSD lengths ##############################
###########################################################################

# Produce figure of predicted TSD lengths.
TE.Names<-c("Ty1", "Ty2", "Ty3", "Ty4")

bed.files1<-list.files(Sys.glob(results.location), full.names=T, pattern="temp_nonredundant.bed$")
bed.files2<-list.files(Sys.glob(results.location), full.names=T, pattern="mapper_nonredundant.bed$")
bed.files3<-list.files(Sys.glob(results.location), full.names=T, pattern="relocate_nonredundant.bed$")
bed.files.rev1<-list.files(Sys.glob(reverse.results.location), full.names=T, pattern="temp_nonredundant.bed$")
bed.files.rev2<-list.files(Sys.glob(reverse.results.location), full.names=T, pattern="mapper_nonredundant.bed$")
bed.files.rev3<-list.files(Sys.glob(reverse.results.location), full.names=T, pattern="relocate_nonredundant.bed$")
all.bed<-c(bed.files1, bed.files.rev1,bed.files2, bed.files.rev2,bed.files3, bed.files.rev3)

# For each bed file
for(i in 1:length(all.bed))
{
	system(paste("grep non- ", all.bed[i], " >> ", output.location, "/testall.bed" ,sep=""))
	dataset <- read.table(paste(output.location, "/testall.bed" ,sep=""), header=FALSE,col.names=c("Chromosome","Start","End", "TE Name", "Score", "Orientation"),sep="\t")
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

new.TEs<-grep("non-reference",dataset$TE.Name)

# Count the number of TEs per sample for each method
newTEs<-paste(dataset$Method[new.TEs],dataset$Sample[new.TEs])
new.TE.counts<-table(newTEs)

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

pdf(paste(output.location, '/singleinsTSDlengths.pdf', sep=""))
par(mfcol=c(4,3), mar=c(2.5,2,1,1), oma=c(4,2.5,2,0))

# ngs_te_mapper TSD lengths
for ( i in 1:length(TE.Names))
{
	TE.tsds<-which(uniq.ngs_te_mapper.insertions$TE == TE.Names[i])
	uniq.ngs_te_mapper.insertions.matrix<-matrix(uniq.ngs_te_mapper.insertions[TE.tsds,9])
	if(TE.Names[i] == "Ty1")
	{
		hist(uniq.ngs_te_mapper.insertions.matrix, xlab="", ylab="Number of insertions", xlim=c(0,15), ylim=c(0,50), main="ngs_te_mapper", breaks=seq(min(uniq.ngs_te_mapper.insertions.matrix)-0.5, max(uniq.ngs_te_mapper.insertions.matrix)+0.5, by=1), col=c(ngs_te_mapper.sr.colour))
	}
	if(TE.Names[i] == "Ty2")
	{
		hist(uniq.ngs_te_mapper.insertions.matrix, xlab="", ylab="Number of insertions", xlim=c(0,15), ylim=c(0,50), main="", breaks=seq(min(uniq.ngs_te_mapper.insertions.matrix)-0.5, max(uniq.ngs_te_mapper.insertions.matrix)+0.5, by=1), col=c(ngs_te_mapper.sr.colour))
	}
	if(TE.Names[i] == "Ty3")
	{
		hist(uniq.ngs_te_mapper.insertions.matrix, xlab="", ylab="Number of insertions", xlim=c(0,15), ylim=c(0,50), main="", breaks=seq(min(uniq.ngs_te_mapper.insertions.matrix)-0.5, max(uniq.ngs_te_mapper.insertions.matrix)+0.5, by=1), col=c(ngs_te_mapper.sr.colour))
	}
	if(TE.Names[i] == "Ty4")
	{
		hist(uniq.ngs_te_mapper.insertions.matrix, xlab="TSD length (bp)", ylab="Number of insertions", xlim=c(0,15), main="", ylim=c(0,50), breaks=seq(min(uniq.ngs_te_mapper.insertions.matrix)-0.5, max(uniq.ngs_te_mapper.insertions.matrix)+0.5, by=1), col=c(ngs_te_mapper.sr.colour))
	}
}

# RelocaTE TSD lengths
for ( i in 1:length(TE.Names))
{
	TE.tsds<-which(uniq.relocate.insertions$TE == TE.Names[i])
	uniq.relocate.insertions.matrix<-matrix(uniq.relocate.insertions[TE.tsds,9])
	if(TE.Names[i] == "Ty1")
	{
		hist(uniq.relocate.insertions.matrix, xlab="",ylab="", xlim=c(0,15), ylim=c(0,50), main="RelocaTE", breaks=seq(min(uniq.relocate.insertions.matrix)-0.5, max(uniq.relocate.insertions.matrix)+0.5, by=1), col=c(relocate.sr.colour))
	}
	if(TE.Names[i] == "Ty2")
	{
		hist(uniq.relocate.insertions.matrix, xlab="",ylab="", xlim=c(0,15), ylim=c(0,50), main="", breaks=seq(min(uniq.relocate.insertions.matrix)-0.5, max(uniq.relocate.insertions.matrix)+0.5, by=1), col=c(relocate.sr.colour))
	}
	if(TE.Names[i] == "Ty3")
	{
		hist(uniq.relocate.insertions.matrix, xlab="",ylab="", xlim=c(0,15), ylim=c(0,50), main="", breaks=seq(min(uniq.relocate.insertions.matrix)-0.5, max(uniq.relocate.insertions.matrix)+0.5, by=1), col=c(relocate.sr.colour))
	}
	if(TE.Names[i] == "Ty4")
	{
		hist(uniq.relocate.insertions.matrix, xlab="TSD length (bp)",ylab="", xlim=c(0,15), ylim=c(0,50), main="", breaks=seq(min(uniq.relocate.insertions.matrix)-0.5, max(uniq.relocate.insertions.matrix)+0.5, by=1), col=c(relocate.sr.colour))
	}
}

# TEMP TSD lengths
for ( i in 1:length(TE.Names))
{
	TE.tsds<-which(uniq.temp.insertions$TE == TE.Names[i])
	uniq.temp.insertions.matrix<-matrix(uniq.temp.insertions[TE.tsds,9])
	if(TE.Names[i] == "Ty1")
	{
		hist(uniq.temp.insertions.matrix, xlab="", ylab="", xlim=c(0,15), ylim=c(0,50),  main="TEMP (SR)", breaks=seq(min(uniq.temp.insertions.matrix)-0.5, max(uniq.temp.insertions.matrix)+0.5, by=1), col=c(temp.sr.colour))
	}
	if(TE.Names[i] == "Ty2")
	{
		hist(matrix(100), xlab="", ylab="", xlim=c(0,15), ylim=c(0,50), col=c(temp.sr.colour), main="")
	}
	if(TE.Names[i] == "Ty3")
	{
		hist(uniq.temp.insertions.matrix, xlab="", ylab="", xlim=c(0,15), ylim=c(0,50), main="", breaks=seq(min(uniq.temp.insertions.matrix)-0.5, max(uniq.temp.insertions.matrix)+0.5, by=1), col=c(temp.sr.colour))
	}
	if(TE.Names[i] == "Ty4")
	{
		hist(uniq.temp.insertions.matrix, xlab="TSD length (bp)", ylab="", xlim=c(0,15), ylim=c(0,50), main="", breaks=seq(min(uniq.temp.insertions.matrix)-0.5, max(uniq.temp.insertions.matrix)+0.5, by=1), col=c(temp.sr.colour))
	}
}

mtext(substitute(paste(italic('Ty4'))), side = 2, line = 0.5, outer = TRUE, at=0.14)
mtext(substitute(paste(italic('Ty3'))), side = 2, line = 0.5, outer = TRUE, at=0.37)
mtext(substitute(paste(italic('Ty2'))), side = 2, line = 0.5, outer = TRUE, at=0.63)
mtext(substitute(paste(italic('Ty1'))), side = 2, line = 0.5, outer = TRUE, at=0.89)

mtext("Predicted TSD lengths (bp)", side = 1, line = 0.5, outer = TRUE)

dev.off()

###########################################################################
############################### Venn diagram ##############################
###########################################################################

# Run external script to calculate "correct" TE predictions
system(paste("bash venn_scoring.sh ", dir, " ", output.location, sep=""))
 
# Read results into R
vennmatrix <- read.delim(paste(output.location, '/vennmatrix.tsv', sep=""), header=FALSE)

# Select for each method (TEMP split by SR and RP)
correctngs <- which(!is.na(vennmatrix[,2]))
correctrelocate <- which(!is.na(vennmatrix[,3]))
correcttempsr <- which(!is.na(vennmatrix[,4]))

correcttemprp <- which(!is.na(vennmatrix[,5]))
correctretroseq <- which(!is.na(vennmatrix[,6]))
correctpopool <- which(!is.na(vennmatrix[,7]))
correcttelocate <- which(!is.na(vennmatrix[,8]))

# Calculate the unique set of insertions detected by each class of method and overall
correctsplitread <- unique(c(correctngs,correctrelocate,correcttempsr))
correctreadpair <- unique(c(correcttemprp, correctretroseq, correctpopool, correcttelocate))
correctall <- unique(c(correctsplitread, correctreadpair))

# Calculate the number of insertions missed by each class of method and overall
missedsplitread <- 598 - length(correctsplitread)
missedreadpair <- 598 - length(correctreadpair)
missedoverall <- 598 - length(correctall)

# Create the diagram for split-read methods
SRvenn <- draw.triple.venn(
area1 = length(correctngs),
area2 = length(correctrelocate),
area3 = length(correcttempsr),
n12 = length(calculate.overlap(x = list(correctngs, correctrelocate))[[3]]),
n13 = length(calculate.overlap(x = list(correctngs, correcttempsr))[[3]]),
n23 = length(calculate.overlap(x = list(correctrelocate, correcttempsr))[[3]]),
n123 = length(calculate.overlap(x = list(correctngs, correctrelocate, correcttempsr))[[1]]),
fill = c(ngs_te_mapper.sr.colour, relocate.sr.colour, temp.sr.colour),
lwd = c(0,0,0),
category = c("ngs_te_mapper", "RelocaTE", "TEMP (Split-read)"),
cat.dist=c(0.1,0.1,0.05),
cat.pos=c(-45,45,180)
)

# Create the diagram for read-pair methods
RPvenn <- draw.quad.venn(
area1 = length(correcttemprp),
area2 = length(correcttelocate),
area3 = length(correctretroseq),
area4 = length(correctpopool),
n12 = length(calculate.overlap(x = list(correcttemprp, correcttelocate))[[3]]),
n13 = length(calculate.overlap(x = list(correcttemprp, correctretroseq))[[3]]),
n14 = length(calculate.overlap(x = list(correcttemprp, correctpopool))[[3]]),
n23 = length(calculate.overlap(x = list(correctretroseq, correcttelocate))[[3]]),
n24 = length(calculate.overlap(x = list(correctpopool, correcttelocate))[[3]]),
n34 = length(calculate.overlap(x = list(correctretroseq, correctpopool))[[3]]),
n123 = length(calculate.overlap(x = list(correcttemprp, correctretroseq, correcttelocate))[[1]]),
n124 = length(calculate.overlap(x = list(correcttemprp, correctpopool, correcttelocate))[[1]]),
n134 = length(calculate.overlap(x = list(correcttemprp, correctretroseq, correctpopool))[[1]]),
n234 = length(calculate.overlap(x = list(correctretroseq, correctpopool, correcttelocate))[[1]]),
n1234 = length(calculate.overlap(x = list(correcttemprp, correctretroseq, correctpopool, correcttelocate))[[1]]),
fill = c(temp.rp.colour, telocate.rp.colour, retroseq.rp.colour, popoolationte.rp.colour),
lwd = c(0,0,0,0),
category = c("TEMP (Read-pair)", "TE-locate", "RetroSeq", "PoPoolationTE"),
cat.dist=c(0.25,0.21,0.1,0.1),
cat.pos=c(-35,20,-10,0)
)

# Create the overall diagram
Overallvenn <- draw.pairwise.venn(
area1 = length(correctsplitread),
area2 = length(correctreadpair),
cross.area = length(calculate.overlap(x = list(correctsplitread, correctreadpair))[[3]]),
fill = c("#ED9D9DCC", "#90ABEDCC"),
lwd = c(0,0),  
category = c("Split-read", "Read-pair"),
euler.d = FALSE,
scaled = FALSE, 
rotation.degree = 180,
cat.dist=0.05
)


# start new page
plot.new()
pdf(paste(output.location, '/venns.pdf', sep=""), height=11, width=7)
# setup layout
gl <- grid.layout(nrow=3, ncol=1)

# setup viewports
vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1)
vp.2 <- viewport(layout.pos.col=1, layout.pos.row=2)
vp.3 <- viewport(layout.pos.col=1, layout.pos.row=3)

# init layout
pushViewport(viewport(layout=gl))
# access the first position
pushViewport(vp.1)

# start new base graphics in first viewport
par(new=TRUE, fig=gridFIG())

grid.text("A.", x = unit(0.15, "npc"), y = unit(0.95, "npc"))
#grid.draw(RPvenn1)
grid.draw(gTree(children=SRvenn, vp=viewport(width=0.5,height=0.95)))
grid.text(missedsplitread, x = unit(0.7, "npc"), y = unit(0.25, "npc"))

# done with the first viewport
popViewport()

# move to the next viewport
pushViewport(vp.2)

grid.text("B.", x = unit(0.15, "npc"), y = unit(0.95, "npc"))
#grid.draw(RPvenn2)
grid.draw(gTree(children=RPvenn, vp=viewport(width=0.55,height=1.1)))
grid.text(missedreadpair, x = unit(0.73, "npc"), y = unit(0.25, "npc"))

popViewport()

# move to the next viewport
pushViewport(vp.3)

grid.text("C.", x = unit(0.15, "npc"), y = unit(0.95, "npc"))
#grid.draw(RPvenn2)
grid.draw(gTree(children=Overallvenn, vp=viewport(width=0.65,height=1.25)))
grid.text(missedoverall, x = unit(0.85, "npc"), y = unit(0.25, "npc"))

# done with this viewport
popViewport(1)
dev.off()

# Remove the empty Rplots.pdf created above
system(paste("rm ", "Rplots.pdf", sep=""))

################## Save workspace and exit the script #####################
save.image(file=paste(output.location, "/workspace.Rdata", sep=""))

###########################################################################