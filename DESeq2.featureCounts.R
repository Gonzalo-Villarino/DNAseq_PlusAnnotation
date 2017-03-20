rm(list=ls())
setwd("~/Documents/NCSU/RNAseq_Mimulus/featureCounts/")
library("DESeq2")
library("edgeR") # u need this for the RPKM calculations
library("RColorBrewer")
library("gplots")

####################################################################################################################################
# Run DESeq2   
####################################################################################################################################

# Import the 6 htseq-count tables to start Diff.Expr.Analysis
# get data
Day00_BR1= read.table(file="0d_BR1_S1.txt",header=FALSE, sep="\t")
colnames(Day00_BR1)[1]<-"gene_short_name"

Day00_BR4= read.table(file="0d_BR4_S10.txt",header=FALSE, sep="\t")
colnames(Day00_BR4)[1]<-"gene_short_name"

Day00_BR5= read.table(file="0d_BR5_S15.txt",header=FALSE, sep="\t")
colnames(Day00_BR5)[1]<-"gene_short_name"

Day00_BR6= read.table(file="0d_BR6_S17.txt",header=FALSE, sep="\t")
colnames(Day00_BR6)[1]<-"gene_short_name"

Day00_BR7= read.table(file="0d_BR7_S19.txt",header=FALSE, sep="\t")
colnames(Day00_BR7)[1]<-"gene_short_name"

Day02d_BR4=read.table(file="2d_BR4_S11.txt",header=FALSE, sep="\t")
colnames(Day02d_BR4)[1]<-"gene_short_name"

Day02d_BR5=read.table(file="2d_BR5_S16.txt",header=FALSE, sep="\t")
colnames(Day02d_BR5)[1]<-"gene_short_name"

Day02d_BR6=read.table(file="2d_BR6_S18.txt",header=FALSE, sep="\t")
colnames(Day02d_BR6)[1]<-"gene_short_name"

Day02d_BR7=read.table(file="2d_BR7_S20.txt",header=FALSE, sep="\t")
colnames(Day02d_BR7)[1]<-"gene_short_name"

Day04_BR2= read.table(file="4d_BR2_S4.txt",header=FALSE, sep="\t")
colnames(Day04_BR2)[1]<-"gene_short_name"

Day04_BR3= read.table(file="4d_BR3_S7.txt",header=FALSE, sep="\t")
colnames(Day04_BR3)[1]<-"gene_short_name"

Day04_BR4= read.table(file="4d_BR4_S12.txt",header=FALSE, sep="\t")
colnames(Day04_BR4)[1]<-"gene_short_name"

Day06_BR1= read.table(file="6d_BR1_S2.txt",header=FALSE, sep="\t")
colnames(Day06_BR1)[1]<-"gene_short_name"

Day06_BR2= read.table(file="6d_BR2_S5.txt",header=FALSE, sep="\t")
colnames(Day06_BR2)[1]<-"gene_short_name"

Day06_BR3= read.table(file="6d_BR3_S8.txt",header=FALSE, sep="\t")
colnames(Day06_BR3)[1]<-"gene_short_name"

Day06_BR4= read.table(file="6d_BR4_S13.txt",header=FALSE, sep="\t")
colnames(Day06_BR4)[1]<-"gene_short_name"

Day07_BR1= read.table(file="7d_BR1_S21.txt",header=FALSE, sep="\t")
colnames(Day07_BR1)[1]<-"gene_short_name"

Day07_BR2= read.table(file="7d_BR2_S22.txt",header=FALSE, sep="\t")
colnames(Day07_BR2)[1]<-"gene_short_name"

Day08_BR1= read.table(file="8d_BR1_S3.txt",header=FALSE, sep="\t")
colnames(Day08_BR1)[1]<-"gene_short_name"

Day08_BR2= read.table(file="8d_BR2_S6.txt",header=FALSE, sep="\t")
colnames(Day08_BR2)[1]<-"gene_short_name"

Day08_BR3= read.table(file="8d_BR3_S9.txt",header=FALSE, sep="\t")
colnames(Day08_BR3)[1]<-"gene_short_name"

Day08_BR4= read.table(file="8d_BR4_S14.txt",header=FALSE, sep="\t")
colnames(Day08_BR4)[1]<-"gene_short_name"




# Merge single htseq-count files into one:
files <- c(
        "0d_BR1_S1.txt",
        "0d_BR4_S10.txt",
        "0d_BR5_S15.txt",
        "0d_BR6_S17.txt",
        "0d_BR7_S19.txt",
        
        "2d_BR4_S11.txt",
        "2d_BR5_S16.txt",
        "2d_BR6_S18.txt",
        "2d_BR7_S20.txt",
        
        "4d_BR2_S4.txt",
        "4d_BR3_S7.txt",
        "4d_BR4_S12.txt",
        
        "6d_BR1_S2.txt",
        "6d_BR2_S5.txt",
        "6d_BR3_S8.txt",
        "6d_BR4_S13.txt",
        
        "7d_BR1_S21.txt",
        "7d_BR2_S22.txt",
        "8d_BR1_S3.txt",
        "8d_BR2_S6.txt",
        "8d_BR3_S9.txt",
        "8d_BR4_S14.txt"
)

Anno=read.table(file="Mguttatus_256_v2.0.annotation_info.txt",header=TRUE,
                sep="\t", comment.char = "", quote = "\"") 
Anno=Anno[,c(-1,-3,-4,-5,-6,-7,-8,-9,-10,-12,-13)]




for (filename in files) {
        AllCountsFinal = read.table(filename,sep="\t")
        AllCountsFinal = merge(AllCountsFinal,Anno,by=c(1))
        #AllCountsFinal = AllCountsFinal[AllCountsFinal[,3]=="protein_coding_gene",]
        #AllCountsFinal = AllCountsFinal[,-3]
        write.table(AllCountsFinal,file=paste(filename,".pro.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
}

#after for loop 3 columns

############### 
#Import files

sampleFiles <- grep("*pro",list.files("~/Documents/NCSU/RNAseq_Mimulus/featureCounts/"),value=TRUE)
sampleFiles

AllCountsFinal=AllCountsFinal[,-3]

#poolA treated
sampleCondition<-c("Day00","Day00","Day00","Day00","Day00","Day02",
                   "Day02","Day02","Day02","Day04","Day04","Day04",
                   "Day06","Day06","Day06","Day06","Day07","Day07",
                   "Day08","Day08","Day08","Day08")
                   
sampleCondition

sampleTable<-data.frame(sampleName=sampleFiles,fileName=sampleFiles, condition=sampleCondition)
sampleTable



ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory="~/Documents/NCSU/RNAseq_Mimulus/featureCounts/",design=~condition)
ddsHTSeq





#Levels in colData are important bc they're used in log calculations; set untreated/control 1st 

colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("Day00","Day02","Day04","Day06",
                                                                          "Day07","Day08"))
ddsHTSeq

#Standard differential expression analysis steps are wrapped into a single function, DESeq.
dds<-DESeq(ddsHTSeq)
str(dds) 
res<-results(dds)
res<-res[order(res$pvalue),] # order results table by the smallest adjusted p value:


#Information about variables and tests were used can be found by calling the function mcols on the results object.
mcols(res)$description
mcols(res,use.names=TRUE)

head(results(dds, addMLE=TRUE),10)

# Results which pass an adjusted p value threshold
resSig <- subset(res, pvalue  < 1) # genes with p-val SMALLER than 0.1
resSig 
dim(resSig)

#write.table(x=resSig, file="BL_YFPs_DESeq2_p-val_0.01")

#####################################################


rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

print(plotPCA(rld, intgroup=c('condition')))
dev.copy(png,'deseq2_pca.png')
dev.off()
