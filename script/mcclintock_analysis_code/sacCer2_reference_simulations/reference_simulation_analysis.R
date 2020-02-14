library(xtable)

# Argument 1 should be the path to the sacCer2_reference_simulations folder
# Argument 2 should be the full path to the desired output location
args<-commandArgs(TRUE)
dir<-args[1]
output.location<-args[2]
dir.create(output.location)

###############################################################
##################### Reference Simulations ###################
###############################################################

tests=array(c("normalreference", "consensusreference", "refcopiesreference", "fullreference", "RMnormalreference", "RMconsensusreference", "RMrefcopiesreference", "RMfullreference"))

nonreftable<-matrix(, nrow = 16, ncol = 10)
reftable<-matrix(, nrow = 16, ncol = 9)

nonreftable[,1]<-c("10X", "10X", "10X", "10X", "10X", "10X", "10X", "10X", "100X", "100X", "100X", "100X", "100X", "100X", "100X", "100X")
nonreftable[,2]<-c("Carr", "Carr", "Carr", "Carr", "RM", "RM", "RM", "RM", "Carr", "Carr", "Carr", "Carr", "RM", "RM", "RM", "RM")
nonreftable[,3]<-c("", "\\tick", "", "\\tick", "", "\\tick", "", "\\tick", "", "\\tick", "", "\\tick", "", "\\tick", "", "\\tick")
nonreftable[,4]<-c("", "", "\\tick", "\\tick", "", "", "\\tick", "\\tick", "", "", "\\tick", "\\tick", "", "", "\\tick", "\\tick")

reftable[,1]<-c("10X", "10X", "10X", "10X", "10X", "10X", "10X", "10X", "100X", "100X", "100X", "100X", "100X", "100X", "100X", "100X")
reftable[,2]<-c("Carr", "Carr", "Carr", "Carr", "RM", "RM", "RM", "RM", "Carr", "Carr", "Carr", "Carr", "RM", "RM", "RM", "RM")
reftable[,3]<-c("", "\\tick", "", "\\tick", "", "\\tick", "", "\\tick", "", "\\tick", "", "\\tick", "", "\\tick", "", "\\tick")
reftable[,4]<-c("", "", "\\tick", "\\tick", "", "", "\\tick", "\\tick", "", "", "\\tick", "\\tick", "", "", "\\tick", "\\tick")

t<-1
for(test in tests)
{
	nonrefmethods=array(c("ngs_te_mapper", "relocate", "temp", "retroseq", "popoolationte", "telocate"))
	refmethods=array(c("ngs_te_mapper", "relocate", "temp", "popoolationte", "telocate"))
	thisnonrefrow10<-vector()
	thisnonrefrow100<-vector()
	for (method in nonrefmethods)
	{
		nonrefmethodcount10<-vector()
		nonrefmethodcount100<-vector()
		for(rep in 1:100)
		{
			nonrefmethodcount10<-cbind(nonrefmethodcount10, system(paste("grep non-reference ", dir, "/", test, "/sacCer2/10X", rep, "simulation/results/10X", rep, "simulation_", method, "_nonredundant.bed | wc | awk '{print $1}'", sep=""), intern=TRUE))
			nonrefmethodcount100<-cbind(nonrefmethodcount100, system(paste("grep non-reference ", dir, "/", test, "/sacCer2/100X", rep, "simulation/results/100X", rep, "simulation_", method, "_nonredundant.bed | wc | awk '{print $1}'", sep=""), intern=TRUE))
		}
		nonrefmethodcount10res<-paste(sprintf("%.2f", round(mean(as.numeric(nonrefmethodcount10)), 2))," \\textpm{} ",sprintf("%.2f", round(sd(nonrefmethodcount10), 2)), sep="")
		nonrefmethodcount100res<-paste(sprintf("%.2f", round(mean(as.numeric(nonrefmethodcount100)), 2))," \\textpm{} ",sprintf("%.2f", round(sd(nonrefmethodcount100), 2)), sep="")
		thisnonrefrow10<-cbind(thisnonrefrow10, nonrefmethodcount10res)
		thisnonrefrow100<-cbind(thisnonrefrow100, nonrefmethodcount100res)
	}
	thisrefrow10<-vector()
	thisrefrow100<-vector()
	for (method in refmethods)
	{
		refmethodcount10<-vector()
		refmethodcount100<-vector()
		for(rep in 1:100)
		{
			refmethodcount10<-cbind(refmethodcount10, system(paste("grep _reference ", dir, "/", test, "/sacCer2/10X", rep, "simulation/results/10X", rep, "simulation_", method, "_nonredundant.bed | wc | awk '{print $1}'", sep=""), intern=TRUE))
			refmethodcount100<-cbind(refmethodcount100, system(paste("grep _reference ", dir, "/", test, "/sacCer2/100X", rep, "simulation/results/100X", rep, "simulation_", method, "_nonredundant.bed | wc | awk '{print $1}'", sep=""), intern=TRUE))
		}
		refmethodcount10res<-paste(sprintf("%.2f", round(mean(as.numeric(refmethodcount10)), 2))," \\textpm{} ",sprintf("%.2f", round(sd(refmethodcount10), 2)), sep="")
		refmethodcount100res<-paste(sprintf("%.2f", round(mean(as.numeric(refmethodcount100)), 2))," \\textpm{} ",sprintf("%.2f", round(sd(refmethodcount100), 2)), sep="")
		thisrefrow10<-cbind(thisrefrow10, refmethodcount10res)
		thisrefrow100<-cbind(thisrefrow100, refmethodcount100res)
	}
	
	nonreftable[t,5:10]<-thisnonrefrow10
	nonreftable[t+8,5:10]<-thisnonrefrow100
	reftable[t,5:9]<-thisrefrow10
	reftable[t+8,5:9]<-thisrefrow100
	t<-t+1
}


latextable1<-xtable(nonreftable)
print(latextable1, file = paste(output.location, '/additional_table1.tex', sep=""), hline.after = c(), include.rownames = FALSE, include.colnames = FALSE, only.contents = TRUE, sanitize.text.function = function(x) x)

latextable1<-xtable(reftable)
print(latextable1, file = paste(output.location, '/additional_table2.tex', sep=""), hline.after = c(), include.rownames = FALSE, include.colnames = FALSE, only.contents = TRUE, sanitize.text.function = function(x) x)

###############################################################