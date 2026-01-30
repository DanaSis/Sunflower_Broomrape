# Rscript Plot.Manhattan.QQ.R trait.assoc.txt trait_name bonferrni_value simplem_value

library(qqman)
library(dplyr) 

options(scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
In=args[1]             
Tag=args[2] 
no_of_SNPs=args[3]
simplem=args[4]          
#kmers=args[5]
print(args)

print("Begin Manhattan plot....")

#import data from GEMMA output folder.
Raw <- read.table(In,header=TRUE) 
#kmers <- read.table(kmers,header=TRUE)


#create new table with only P-values from wald test.
Data <- cbind(Raw[,1:3],Raw$p_wald) %>%
  mutate(logged_p = -log10(Raw$p_wald))  # Calculate the logged p-value too

max_y_limit <- max(Data$logged_p)


#kmers <- kmers[c(1,4,2,3)]
Data <- Data[c(1,2,3,4)]
head(Data)
#head(kmers)
colnames(Data) <- c("CHR","SNP","BP","P")
#colnames(kmers) <- c("CHR","SNP","BP","P")
#Datak <- rbind(Data, kmers) 
#change colunm names.
colnames(Data) <- c("CHR","SNP","BP","P")
Data$CHR <- as.numeric(gsub("Ha412HOChr", "", Data$CHR))
#Data$CHR <- as.numeric(gsub("HanXRQChr", "", Data$CHR))
#Data$CHR <- as.numeric(gsub("Ha", "", Data$CHR))


#remove all letters in this columns and change class of the chromosome number to character
#Data$CHR <- as.numeric(str_replace_all(as.character(Data$CHR), "[a-z,A-Z]",""))



#SNPs=c("chr4H_554293917","chr7H_85388149","chr6H_373658496","chr5H_538244017","chr4H_518117317","chr3H_17268612","chr3H_17298634","chr3H_17341304","chr1H_506121081")


#define name for manhattan plot
Man=paste("ManP_",Tag,".jpg",sep="")


#create threshold bonferoni line
no_of_SNPs = as.numeric(no_of_SNPs)
simplem = as.numeric(simplem)

#sug=-log10(0.05/716824) # from the Inferred Meff of simplem output
#sug = -log10(1e-5) #suggestive for GWAS in general
bon =-log10(0.05/no_of_SNPs) # 0.05/no. of SNPs
sim = -log10(simplem) # simplem
#highlights = kmers [[2]]

#create manhattan plot (NOTE:	make sure there are no NAs in the CHR column -no extra unidentified chromosomes before plotting)

#man=paste(Tag,"_ManP.pdf",sep="")
jpeg(paste("ManP_",Tag,".jpeg",sep=""),2000, height = 1000, units = "px", quality = 100)
manhattan(Data, main=Tag,
                 logp = T,
                 ylab = paste("-log10(p)"),
				 ylim=c(0,15),
				 #ylim = c(0,max_y_limit),
                 #suggestiveline =p.adjust(data[,4],method="fdr"), 
                 cex = 2.5, size = 0.05, 
		 cex.axis = 2,
		 cex.main=4,
		 cex.lab=1.5,
			#col = c("plum3","lightgoldenrod"), #snp gwas
			col = c("cadetblue","bisque4"), #kmer gwas
		 #highlight = highlights,
		 genomewideline = sim,
		 suggestiveline = bon)
		 #chrlabs=as.character(c(1:7))) 

#man=paste(Tag,"_ManP.pdf",sep="")
#pdf(filename=man, width=10,height=7)
#manhattan(Data,
 #         main=Tag,
  #        col=c("mediumpurple4","gold2"),
   #       suggestiveline=-log10(7.702012e-08),
    #      genomewideline=-log10(2.147942e-07),
     #     annotatePval = 0.01)


#this is importand to export the plot to the directory- if not written, file size will be 0
dev.off()

print("Manhattan plot - DONE !")

print("Begin plotting Q-Q Plot...") 
#create Q-Q plot

QQ=paste("QQ_",Tag,".jpeg",sep="")
#1.prepare jpg file
jpeg(filename =QQ, width = 500, height = 500, units = "px", quality = 100)

#2.create plot
qq(Data$P)

#3.export
dev.off()

print("Q-Q Plot - DONE !")
