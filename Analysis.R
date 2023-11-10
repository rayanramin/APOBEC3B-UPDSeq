#### This is the main code for the analysis and plotting of the data

## Samples.xlsx lists of all of the used Samples and their accession numbers

## 4 A3A samples are already published 
# A3A_F6 SRR6924522	PRJNA448166 
# A3A_H1 SRR9864913	PRJNA448166
# A3A_A SRR17822878	PRJNA801888
# A3A_G SRR17822877	PRJNA801888

## 3 A3B-CTD + 3 Empty vector control samples
## 6 A3B-full + 6 Empty vector control samples

## ==> run the alignment and bam-readcount for the samples
# bam-readcount -w0 -f <ref.fa> <bam_file> | \
# awk -F ":|\t|=" 'BEGIN {OFS = "\t"}; {print $1, $2, $3 , $4, $21 , $35, $49 , $63}' > out.txt

############################################
##########      NDC2 Analysis     ##########
############################################

# ## Sequence Alignment and extraction of depth of coverage

# # Load necessary modules
# module load samtools/1.9 bedtools/2.25.0 r/3.2.3 bwa/0.7.17 bamtools/2.5.1
# # Index reference DNA files
# bwa index -p gen_ref ./BH214V4.fa
# bwa index -p plas_ref ./plasmid_reference.fa
# # Align and exclude vector-aligned reads
# bwa mem -t 4 plas_ref ./data/fastq/A3B-I-3/A3BS3R1.fastq.gz ./data/fastq/A3B-I-3/A3BS3R2.fastq.gz | samtools view -b -f 12 -o plas_excl_align_A3BI3.bam
# samtools sort -n -l 6 -@ 4 plas_excl_align_A3BI3.bam -o plas_excl_sort_A3BI3.bam
# # Convert BAM to FASTQ 
# bedtools bamtofastq -i plas_excl_sort_A3BI3.bam -fq plas_excl_A3BI3_R1.fq -fq2 plas_excl_A3BI3_R2.fq
# # Align reads to genomic DNA
# bwa mem -t 4 gen_ref plas_excl_A3BI3_R1.fq plas_excl_A3BI3_R2.fq | samtools view -b -o A3BI3.aligned.bam
# samtools sort -l 6 -@ 4 A3BI3.aligned.bam -o A3BI3.sort.bam
# # Extract depth of coverage
## note, it's important to use the -aa option
# samtools depth -aa -m 100000 A3BI3.sort.bam > A3BI3.depth

### Example: A3BI3 vs Negative control EVI3

### Import depth of coverage files to R
library(RamPack)
indir <- "/Path/to/depth/files/"

list.files(indir,pattern=".depth",full.names = T)
# "A3BI3","EVI3"

A3BI3_DEPTH <- read.table(paste0(indir,"A3BI3.depth"),col.names=c("Chr", "POS", "COV"))
EVI3_DEPTH <- read.table(paste0(indir,"EVI3.depth"),col.names=c("Chr", "POS", "COV"))

# Calculate NDC2 and five sigma standard deviation (FSD) for data and plot it
w = 1e5 ; b = 1e4 ; mavn = 120;
A3BI3_DEPTH$NDC <- RamPack::NDC(x = A3BI3_DEPTH$COV, y = EVI3_DEPTH$COV, w = w, b = b, mavn = mavn,NDC_Version=2)

FSD <- sd(A3BI3_DEPTH$NDC, na.rm=TRUE)*5

### plot NDC
options(bitmapType='cairo')
png(filename = paste0(outdir,"A3BI3_NDC2.png") ,height = 6, width = 24, units = "in", res = 400)
par( mai = c(1,1,0.7,0.2))
xl <- pretty(A3BI3_DEPTH$POS, n = 10)/10^6
plot(A3BI3_DEPTH$POS/1e6, A3BI3_DEPTH$NDC, type = "l",main = "NDC2: A3BI3 vs EVI3" ,axes=F, cex.main =2.5 , xlab = "", ylab="")
axis(2,cex.axis=2)
axis (1, cex.lab=2 ,cex.axis = 2,cex=2 , at =xl , labels = paste(xl, "Mbp", sep = ""), las=1)
abline(h= FSD , col="red")
abline(h= -FSD , col="red")
title( ylab = "Relative Coverage", xlab ="Genomic Position", cex.lab=2)
box()
dev.off()

# Find regions of significant enrichment 
T_SPD <- R.utils::seqToIntervals(A3BI3_DEPTH$POS[A3BI3_DEPTH$NDC > FSD]) %>% as.data.frame
#########
#Export list of uracil peaks into bed file 
fixSPD = function(x, chrName){                           
z <- x
z[,3] <-z[,2]
z[,2]<-z[,1]
z[,1] <- chrName
names(z) <- c("chr.", "from","to")
return(z)
}
#########
bedTools.2in = function(functionstring="bedIntersect",bed1,bed2,opt.string=""){
 #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=read.table(out,header=F, sep="\t")
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}
#########
fixSPD(x=T_SPD,"BH214") -> T_SPD
## save the peaks as bed file
write.table(T_SPD, paste0(topdir,"A3BI3_NDC2",".peaks.bed"), quote=F, sep="\t", row.names=F, col.names=F)
#########
# Annotation 
annotations <- read.table(paste0(topdir,"BH214.gff"), sep="\t") #import list of annotations
## bedtools is required for this step
j <- bedTools.2in("bedtools intersect", annotations,T_SPD)
## save the list of overlapping genes
write.table(j, paste(topdir, "A3BI3_NDC2",".genes.txt", sep="_"), quote=F, sep="\t", row.names=F, col.names=F)
############################################
## Run the above code on all samples to generate the list of uracil peaks
## then run the following code to extract the intersection of peaks from samples
## using bedtools in command line:

#$ bedtools intersect -a A3Bctd_rep1.bed -b A3Bctd_rep2.bed A3Bctd_rep3.bed > intersect_345A3B-CTD.bed
#$ bedtools intersect -a A3Bfull_rep1.bed -b A3Bfull_rep2.bed A3Bfull_rep3.bed > intersect123A3Bfull.bed
#$ bedtools intersect -a A3A_rep1.bed -b A3A_rep2.bed A3A_rep3.bed A3A_rep4.bed > intersectA1234.bed

##==> The bedfiles are provided in the folder "Uracilation_Peaks_bed"

### Figure 1A (Barcode plot for A3A, A3B-CTD and A3B-full intersect uracil peaks)

indir <- "/path/to/bed/files/"
# import the intersection of uracilation peaks
# common peaks in A3B_full
A3B_FULL <- read.delim(paste0(indir,"intersect123A3Bfull.bed"), header = FALSE)
# common peaks in A3A
A3A <- read.delim(paste0(indir,"intersectA1234.bed"), header = FALSE)
# common peaks in A3B-CTD
A3B_CTD <- read.delim(paste0(indir,"intersect_345A3B-CTD.bed"), header = FALSE)
# common peaks between A3A and A3B-full (asterisk in the figure)
intersect1 <- read.delim(paste0(indir,"intersect_A3A4rep_A3B-FL.bed"), header = FALSE)
# common peaks between A3A and A3B-CTD (asterisk in the figure)
intersect2 <- read.delim(paste0(indir,"intersect_A3A_A3B-CTD_A3B-FL.bed"), header = FALSE)

j = 3000 # to control the thickness of the barcode lines
ggplot(A3B_FULL) +
geom_rect(aes(xmin = (V2 -j)/1e6 , xmax = (V3 +j)/1e6 , ymin = 0.05 , ymax = 0.75 ),fill="black") +
geom_rect(A3A ,mapping=aes(xmin = (V2 -j)/1e6 , xmax = (V3 +j)/1e6 , ymin = 1 , ymax = 1.75),fill="black") +
geom_rect(A3B_CTD ,mapping=aes( xmin =(V2-j)/1e6 , xmax = (V2+j)/1e6 , ymin = 2 , ymax = 2.75 ),fill="black") +
annotate("text" , x = -.2 , y = 0.5 , label = "A3B-FULL" , angle = 90)+
annotate("text" , x = -.2 , y = 1.5 , label = "A3A" , angle = 90)+
annotate("text" , x = -.2 , y = 2.4 , label = "A3B-CTD" , angle = 90)+
labs(title = "All intersect peaks of A3B-CTD, A3A, and A3B-FULL" ,x = "Genomic Position (Mbp)")+
annotate("text" , x = c(rowMeans(intersect1[,c(2:3)]))/1e6 , y = 0.8 , label = "*" , size  =6, color= "#4FC978" )+
annotate("text" , x = c(rowMeans(intersect2[,c(2:3)]))/1e6 , y = 2.75 , label = "*" , size  =6, color= "#c941f2" )+
scale_x_continuous(breaks = 0:4)+
coord_cartesian(ylim = c(0,2.8) , xlim = c(-.2,4.7)) + 
theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
	plot.margin = unit(c(0.1,.2,.4,.1) ,units = "in"),
	axis.text.x   = element_text(size=14),
	axis.ticks.x.bottom = element_line(size =1),
	axis.ticks.length =  unit(c(.05) ,"in") ,
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y =element_blank(),
	axis.title.y = element_blank(),
	axis.title.x = element_text(size = 15),
	plot.background = element_blank(),
	panel.background = element_blank(), 
	legend.position = "") -> g.1a
ggsave(filename= paste0(outdir,"Fig1A_Barcode.pdf"),plot=g.1a, height = 4, width = 8, units = "in" )




############################################
# For the Uracilation Index analysis 3 data.frames are needed:
# 1. Samples: a data.frame listing the samples
# 2. X: a data.frame containing the uracilation fractions for all positions in all samples
# 3. POSData: result of survey-hairpin on BH214 genome #  Result out of ApoHP (https://github.com/alangenb/ApoHP) 
############################################
# 1- list of  samples
Samples <- data.frame(
	stype=c(rep("A3A",4),rep("A3B_CTD",3),rep("A3B_full",6),rep("EV",9)),
	samp_rep=c(paste0("A3A_",1:4),paste0("A3B_CTD_",1:3),paste0("A3B_full_",1:6),paste0("EV_",1:9)),
	samp=c("A3A_F6","A3A_H1","A3A_A","A3A_G","A3BI3","A3BI4","A3BI5",
	"A3Bfl_B3S1","A3Bfl_B3S2","A3Bfl_B3S3","A3Bfl_R1_F8S6","A3Bfl_R1_F10S4","A3Bfl_R1_F12S3",
	"EV1_S4","EV3_S5","EV3S2","EV4S1","EV5_S6","EV5S5","EVI3","EVI4","EVI5"))
####  
# 2- X (uracilation fractions calculated from readcounts)
# example bam-readcount: 
# bam-readcount -w0 -f <ref.fa> <bam_file> | awk -F ":|\t|=" 'BEGIN {OFS = "\t"}; {print $1, $2, $3 , $4, $21 , $35, $49 , $63}' > out.txt

topdir <- "/path/to/your/data/"
indir = "/path/to/readcount/files/" # this is where you have your readcount files
files <- list.files(indir,pattern=".readcount.txt",full.names = T)
##
all_readcounts = list()
for(i in 1:length(files)){
	print(paste("loading ",Samples$samp_rep[i]));
	infile = files[grep(Samples$samp[i],files)];
	RamPack::read.readcount(infile) -> 	all_readcounts[[Samples$samp[i]]]
}
####
# function to calculate uracilation fraction
# also available from the RamPack R package
# devtools::install_github("rayanramin/RamPack")
UX <- function(Y) {
    X <- Y
    # defining conditions    
    cond1 <- logical(nrow(X)) ; cond1 <- X$REF == "C"
    cond2 <- logical(nrow(X)) ; cond2 <- X$REF == "G"
    # C to T changes only
    t <- ifelse(cond1, as.integer(X$T), 0) / X$COV
    # G to A changes only
    a <- ifelse(cond2, as.integer(X$A), 0) / X$COV
    # All changes
    c(a+t) 
}
####
# get the uracilation fractions for all samples
library(dplyr); library(RamPack)
X = NULL
lapply(1:length(all_readcounts) ,
	function(x){
		all_readcounts[[x]] %>% 
		filter(REF %in% c("C","G"))	%>%
		mutate(samp =Samples$samp[x],
		U = UX(.)) %>%
		select(POS,U,samp)
	}
) %>%  bind_rows() -> X
##############################################
# 3- POSData
# import the result of ApoHP survey_hairpins on BH214 genome
POSData <- read.table(paste0(topdir, "AllC_BH214.txt" ), header = T , sep = "\t")
names(POSData)[1] <- "POS"
## stem strength bins
POSData$ssGroup <-  cut(POSData$ss , breaks=c(-1,7,9,11,13,15,17,19,Inf), labels=c("0-7","8-9","10-11","12-13","14-15","16-17","18-19", "20+"))
######## Add LDST/LGST info to POSData
POSData %>% {ifelse( .$ref == "2" , 1 , -1 )} -> POSData$strand
POSData %>% {ifelse((.$POS > 1664141 & .$POS < 3906549 ), "Left", "Right" )} %>% as.factor() -> POSData$Replication
POSData %>% {ifelse(((.$Replication == "Right" & .$strand == 1) | (.$Replication == "Left" & .$strand == -1)), "LGST" , "LDST")} %>% as.factor() -> POSData$LDLGST
########
### save the data
save(Samples,X,POSData,file=paste0(topdir,"Uracilation_Data.RData"))
##############################################
##############################################
## Perform the analysis
load(paste0(topdir,"Uracilation_Data.RData"))
library(dplyr); library(stringr); library(tidyr);
library(RamPack);library(ggplot2);library(Rmisc);
library(magrittr); library(ggrepel)
############
pval2asterisk <- function(p){
  cut(p, breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), 
                labels = c("****", "***", "**","*", "n.s."), right = FALSE)
}
#############
# Load the  BH214 genome, necessary for Fig 2C, 2E,2F and Fig3
# install.packages("BSgenome.Ecoli.WayneState.BH214V4_1.0.0.tar.gz",repos=NULL,type="source")

ref_genome <- "BSgenome.Ecoli.WayneState.BH214V4"
library(ref_genome, character.only = TRUE)
#############
# function to get the sequence of a hairpin
get_hairpin_seq <- function(D,REF=BH214$BH214V4, stem_length=1, onlyLoop=FALSE){
	# D is the POSData data.frame and requires the following columns:
	# POS, looplen, looppos, ref 
	D %>% mutate(roc = looplen-looppos, loc = looppos-1) %>%
      mutate(st = ifelse(ref==2, (POS-loc), (POS-roc)) , en = ifelse(ref==2, POS+roc, POS+loc)) %>%
      mutate(stem_st = st - stem_length , stem_en = en + stem_length ) -> D
	D$seq <- NA

	## to deal with circular genome
	# beginning
	D$seq[which(D$stem_st <=0)] <- lapply(which(D$stem_st <=0), function(i) paste0(REF[(length(REF)+D$stem_st[i]):length(REF)] , REF[1:D$stem_en[i]])) %>% unlist %>% toupper
	## end
	D$seq[which(D$stem_en >length(REF))] <- lapply(which(D$stem_en >length(REF)), function(i) paste0(REF[(D$stem_st[i]):length(REF)] , REF[1:(D$stem_en[i]-length(REF))])) %>% unlist %>% toupper
	## middle of the sequence
	D$seq[ which(D$stem_st >0 & D$stem_en <=length(REF))] <- lapply(which(D$stem_st >0 & D$stem_en <=length(REF)), function(i) paste0(REF[(D$stem_st[i]):D$stem_en[i]])) %>% unlist %>% toupper
	## Reverse complement when the reference is G
	D$seq[D$ref==3] <- lapply(which(D$ref==3), function(i) seqinr::c2s(rev(seqinr::comp(seqinr::s2c(D$seq[i])))) ) %>% unlist %>% toupper

	### Sanity check
	if (any(!(substr(D$seq,stem_length+D$looppos,stem_length+D$looppos)=="C"))) {
	stop("something went wrong")
	} else {
		SEQS <-  lapply(D$seq, seqinr::s2c) 
		lapply(1:length(SEQS), function(i) {
			# denote the targeted C
			SEQS[[i]][stem_length+D$looppos[i]] <- ".C." ; 
			# add parenthesis to denote the loop
    		SEQS[[i]][1+stem_length] <-  paste0('(',SEQS[[i]][1+stem_length] ) ; 
    		SEQS[[i]][stem_length+D$looplen[i]] <-  paste0(SEQS[[i]][stem_length+D$looplen[i]],')' );
    		seqinr::c2s(SEQS[[i]]) }) %>%  unlist -> out
	}
	##		
	if(onlyLoop){
	stringr::str_extract_all(out,"\\(.*?(.C.)?\\)|\\.C\\.",simplify = TRUE) %>%
	apply(. ,1,paste,collapse="") -> out
	}
	return(out)
}



#############
# Analysis and plotting
##############

####################################################################################
### Figure 1B
## normalized replication bias

# per sample
left_join(filter(POSData, minus0==4),X,by="POS") %>% 
left_join(.,Samples, by="samp") %>% 
summarySE(measurevar="U", groupvars=c("stype","samp","LDLGST")) %>%  
mutate(UI=U*1e3) %>%
# normalized per sample  and measure over sample type 
select(., -sd,-se,-ci,-N) %>% 
mutate(UN = ave(UI, samp, FUN=function(x)  x/mean(x))) -> tmp1b 
summarySE(tmp1b,measurevar="UN", groupvars=c("stype","LDLGST")) -> Pdata.1b
# add labels for the plot
Pdata.1b %$% ifelse(LDLGST=="LDST",NA,ave(UN, stype, FUN=function(x)  paste0(round((x/min(x)),2),"X")) ) -> Pdata.1b$lab

ggplot(Pdata.1b,aes(x = stype , y = UN ,fill = LDLGST  )) + 
  geom_hline(yintercept = 1 , lty = "dotted" , lwd = 1.5) +
  geom_bar( stat = "identity" , position = position_dodge2(width = 1 ,padding = 0 )  , lwd =.5, col="darkgrey" )+
  geom_errorbar(inherit.aes = F , mapping = aes(x = stype , ymin=UN-sd, ymax=UN+sd , group = LDLGST) , col = "black" , width = 0.2, size=1 , position = position_dodge(.9 )) + # this is messy but it works!
  theme_classic() + ylab("Normalized Uracilation Index" ) + 
  labs(x="", title =  "Replication Bias at TC context" , fill = ""  )  + 
  geom_text(aes(y=UN-.1,label=lab), position=position_dodge(.9),
  size=7, color="white",fontface="bold")+
  scale_y_continuous(breaks = seq(0,1.4,.2)) + 
  coord_cartesian(ylim = c(.4, 1.7)) + 
  scale_x_discrete(labels = c("A3A" , "A3B-CTD" , "A3B-full","EV")) + 
  scale_fill_manual(breaks =c("LDST","LGST"), 
  values = c("#a0d1fd","#6888a3")) +
  theme(axis.text = element_text(size = 23 , colour = "Black") , 
  axis.title.y = element_text(size =26) , 
  axis.text.y = element_text(size = 22) , 
  title = element_text(size =24) , 
  axis.text.x = element_text(face = "bold", vjust = .6) , 
  legend.text =  element_text(size = 25) , 
  legend.direction = "horizontal",
  legend.position = c(0.25,0.94)) +
  guides(fill = guide_legend(keywidth=unit(.7,"in"), keyheight =unit(.7,"in"))) -> g.1b
#############
## p value calculation:
#############
filter(tmp1b,stype=="A3A") %>%
t.test(UN ~ LDLGST, data = ., paired = TRUE)  %$% pval2asterisk(p.value) -> pval_A3A
#############
filter(tmp1b,stype=="A3B_CTD") %>%
t.test(UN ~ LDLGST, data = ., paired = TRUE)  %$% pval2asterisk(p.value) -> pval_A3Bctd
#############
filter(tmp1b,stype=="A3B_full") %>%
t.test(UN ~ LDLGST, data = ., paired = TRUE) %$% pval2asterisk(p.value) -> pval_A3Bfull
#############
g.1b + annotate("text", x = "A3A", y = 1.35, label = pval_A3A, size = 8) + 
geom_segment(aes(x = 0.7, xend = 1.3, y = 1.33, yend = 1.33), color = "black", size = 2) +
annotate("text", x = "A3B_CTD", y = 1.44, label = pval_A3Bctd, size = 8) +
geom_segment(aes(x = 1.7, xend = 2.3, y = 1.42, yend = 1.42), color = "black", size = 2) +
annotate("text", x = "A3B_full", y = 1.47, label = pval_A3Bfull, size = 8) +
geom_segment(aes(x = 2.7, xend = 3.3, y = 1.45, yend = 1.45), color = "black", size = 2) -> g.1b_P

outdir <- paste0(topdir,"Figures/");dir.create(outdir, showWarnings = FALSE)
ggsave(filename= paste0(outdir,"Fig1B_RepBias.pdf"),plot=g.1b_P, height = 10, width = 10, units = "in" )

####################################################################################
####################################################################################
### Figure 2A
## Normalized Uracilation in strong short hairpins
POSData %<>% mutate(hairpin =  {(ss >= 15 & looplen %in% 3:5 & looppos>0 & looppos <= looplen)})  

### extract the data per sample
filter(X , U < 0.8) %>%  
left_join( select( filter(POSData,minus0==4) , POS,  hairpin) , .) %>%  
left_join(.,Samples, by="samp") %>%
# normalized per sample  and measure over sample type 
summarySE(.,measurevar="U", groupvars=c( "stype","samp","hairpin")) %>% 
mutate(UI=U*1e3) %>% 
# normalize so the non-hairpin is 1
mutate(UN = ave(UI, samp, FUN=function(x)  x/x[1])) -> tmp2a
summarySE(tmp2a,  measurevar="UN", groupvars=c( "stype","hairpin")) -> Pdata.2a
######
Pdata.2a %$% ifelse(!hairpin | stype=="EV",NA ,ave(UN, stype, FUN=function(x)  paste0(round((x/min(x)),1),"X")) )  -> Pdata.2a$lab

ggplot(Pdata.2a,aes(x = stype , y = UN ,fill = hairpin  )) + 
  geom_hline(yintercept = 1 , lty = 11 , lwd = 1) +
  geom_bar( stat = "identity" , position = position_dodge2(width = 1 ,padding = 0 )  , lwd =0.5, color="black" )+
  geom_errorbar(inherit.aes = F , mapping = aes(x = stype , ymin=UN-sd, ymax=UN+sd , group = hairpin) , col = "black" , width = 0.2, size=1 , position = position_dodge(.9 )) + # this is messy but it works!
  theme_classic() + labs(y="Normalized Uracilation Index" ) + 
  labs(x="", title =  "Uracilation at Strong Short Hairpins" , fill = ""  )  + 
  geom_text(aes(y=UN-sd-1 ,label=lab), position=position_dodge(.9),
  size=7, color="white",fontface="bold")+
  scale_y_continuous(expand=c(0,0), breaks = seq(0,11,2)) + 
  coord_cartesian(ylim = c(0, 13.5)) + 
  scale_x_discrete(labels = c("A3A" , "A3B-CTD" , "A3B-full","EV")) + 
  scale_fill_manual(breaks =c(T,F),labels= c("Strong Short Hairpins","Rest of the Genome"),
  values = c("#9f001d","#c79278")) +
  theme(axis.text = element_text(size = 23 , colour = "Black") , 
  axis.title.y = element_text(size =26) , 
  axis.text.y = element_text(size = 22) , 
  title = element_text(size =24) , 
  axis.text.x = element_text(face = "bold", vjust = .6) , 
  legend.text =  element_text(size = 25) , 
  legend.direction = "vertical",
  legend.position = c(0.65,0.85)) +
  guides(fill = guide_legend(override.aes = list(linetype = 0),keywidth=unit(.7,"in"), keyheight =unit(.7,"in"))) -> g.2a

#############
# p value calculation:
#############
filter(tmp2a,stype=="A3A") %>%
t.test(UN ~ hairpin, data = ., paired = TRUE)  %$% pval2asterisk(p.value) -> pval_A3A
#############
filter(tmp2a,stype=="A3B_CTD") %>%
t.test(UN ~ hairpin, data = ., paired = TRUE)  %$% pval2asterisk(p.value) -> pval_A3Bctd
#############
filter(tmp2a,stype=="A3B_full") %>%
t.test(UN ~ hairpin, data = ., paired = TRUE) %$% pval2asterisk(p.value) -> pval_A3Bfull
#############
g.2a + annotate("text", x = "A3A", y = 12.66, label = pval_A3A , size = 8) + 
geom_segment(aes(x = 0.7, xend = 1.3, y = 12.56, yend = 12.56), color = "black", size = 2) +
annotate("text", x = "A3B_CTD", y = 7.68, label = pval_A3Bctd , size = 8) +
geom_segment(aes(x = 1.7, xend = 2.3, y = 7.58, yend = 7.58), color = "black", size = 2) +
annotate("text", x = "A3B_full", y = 7.11, label = pval_A3Bfull, size = 8) +
geom_segment(aes(x = 2.7, xend = 3.3, y = 7.01, yend = 7.01), color = "black", size = 2) -> g.2a_P
ggsave(filename= paste0(outdir,"Fig2A_StrongHP.pdf"),plot=g.2a_P, height = 10, width = 10, units = "in" ) 

####################################################################################
### Figure 2B
## Uracilation by Stem Strength

filter(X , U < 0.8) %>%  
left_join( select( filter(POSData,minus0==4) , POS,  ssGroup) , .) %>%  
left_join(.,Samples, by="samp") %>%
### calculate the stats
summarySE(. , measurevar="U", groupvars=c("stype","samp","ssGroup")) %>% 
mutate(U=U*1e3) %>% 
mutate(UN = ave(U, samp, FUN=function(x)  x/x[1])) %>% # normalize so the ss-Group 0-7
summarySE(  measurevar="UN", groupvars=c( "stype","ssGroup")) -> Pdata.2b
######
mutate(Pdata.2b,stype= factor(stype, levels=c("A3A","A3B_CTD", "A3B_full","EV"))) %>% 
ggplot() +
geom_hline(yintercept = 1 , lty = 11 , lwd = 1) +
geom_bar(aes(x=stype, y=UN, fill=ssGroup),alpha=1,lwd = 0.25 , col="black", stat= "identity", position=position_dodge())+
theme_classic()+
scale_fill_brewer(palette = "YlOrRd") + 
labs(title="Uracilation by Stem Strength", x="",y="Normalized Uracilation Index" , fill= "Hairpin\nStem \nStrength")+
scale_y_continuous(expand = c(0,0), breaks=seq(1,6,1))+
geom_errorbar(inherit.aes = F , mapping = aes(x = stype ,group = ssGroup,ymin=UN-sd, ymax=UN+sd), width=.25,linewidth=.5, position=position_dodge(.9)) +
theme( plot.title = element_text(size = 25), 
axis.text.x = element_text(face = "bold",color="black",size = 24), 
axis.text.y  = element_text(size = 25,color="black"), 
axis.title  = element_text(size=28) , 
legend.title = element_text(size = 25), 
legend.text = element_text(size=24) , 
legend.position =   c(.85, .72)) +
coord_cartesian(ylim=c(0,7)) -> g.2b
# ggsave(filename= paste0(outdir,"Fig2B_HairpinSSGroup.pdf"),plot=g.2b, height = 10, width = 12, units = "in" )
#################
# p value calculation:
#################
filter(X , U < 0.8) %>%  
left_join( select( filter(POSData,minus0==4) , POS,  ssGroup) , .) %>%  
left_join(.,Samples, by="samp") %>%
### calculate the stats
summarySE(. , measurevar="U", groupvars=c("stype","samp","ssGroup")) %>% 
mutate(U=U*1e3) %>% 
mutate(UN = ave(U, samp, FUN=function(x)  x/x[1])) -> Data.2b
# A3A: 0-7 vs 20+
filter(Data.2b,stype=="A3A",ssGroup %in% c("0-7","20+")) %>%
t.test(UN ~ ssGroup, data = ., paired = TRUE)  %$% pval2asterisk(p.value) -> pval_A3A
# A3B_CTD: 0-7 vs 20+
filter(Data.2b,stype=="A3B_CTD",ssGroup %in% c("0-7","20+")) %>%
t.test(UN ~ ssGroup, data = ., paired = TRUE)  %$% pval2asterisk(p.value) -> pval_A3B_ctd
# A3B_full: 0-7 vs 20+
filter(Data.2b,stype=="A3B_full",ssGroup %in% c("0-7","20+")) %>%
t.test(UN ~ ssGroup, data = ., paired = TRUE)  %$% pval2asterisk(p.value) -> pval_A3B_full
####
g.2b +
geom_segment(aes(x = 0.6, xend = 1.4, y = 7.2, yend = 7.2), color = "black", size = 1) +
geom_segment(aes(x = 0.6, xend = 0.6, y = 1.3, yend = 7.2), color = "black", size = 1) +
geom_segment(aes(x = 1.4, xend = 1.4, y = 7, yend = 7.2), color = "black", size = 1) +
annotate("text", x = "A3A", y = 7.3, label = pval_A3A, size = 10)+
geom_segment(aes(x = 1.6, xend = 2.4, y = 4.2, yend = 4.2), color = "black", size = 1) +
geom_segment(aes(x = 1.6, xend = 1.6, y = 1.3, yend = 4.2), color = "black", size = 1) +
geom_segment(aes(x = 2.4, xend = 2.4, y = 4, yend = 4.2), color = "black", size = 1) +
annotate("text", x = "A3B_CTD", y = 4.3, label = pval_A3B_ctd, size = 10) +
geom_segment(aes(x = 2.6, xend = 3.4, y = 4.2, yend = 4.2), color = "black", size = 1) +
geom_segment(aes(x = 2.6, xend = 2.6, y = 1.3, yend = 4.2), color = "black", size = 1) +
geom_segment(aes(x = 3.4, xend = 3.4, y = 4, yend = 4.2), color = "black", size = 1) +
annotate("text", x = "A3B_full", y = 4.3, label = pval_A3B_full, size = 10) -> g.2b_P
ggsave(filename= paste0(outdir,"Fig2B_HairpinSSGroup_pvals.pdf"),plot=g.2b_P, height = 10, width = 12, units = "in" )

####################################################################################
### Figure 2C
## correlation among samples based on hairpin targeting

## Filter POSData for the hairpin positions
filter(POSData , looplen %in% 3:6 & looppos >0 & looppos <= looplen & ss >=10) -> POSX
## extract the sequence of the hairpins
POSX %<>% mutate(loopseq = get_hairpin_seq(.,REF=BH214$BH214V4,onlyLoop=T)) 

# count rows by group
group_by(POSX,looplen,looppos,loopseq) %>% 
dplyr::summarise(n_total = n()) -> SX 

select(POSX, POS, looplen, looppos, loopseq) %>% 
left_join(filter(X , U < 0.8))  %>% 
left_join(.,Samples, by="samp") -> POSX_U

## calculate the uracilation fraction per sample
filter(POSX_U, stype != "EV") %>%
mutate(POSX_U, UI = U*1e3) %>%
summarySE(., measurevar="UI", groupvars=c("samp_rep","loopseq")) %>%
left_join(.,SX, by="loopseq") -> SX_U_sample

## calculate correlations
SX_U_sample %>%  
filter(n_total >= 10) %>%
select(samp_rep, loopseq, UI) %>% 
spread(key="samp_rep", value="UI") %>% 
select(-1) %>%  cor( use = "complete.obs") -> SX_U_sample_cor

#####
reshape2::melt(SX_U_sample_cor) -> Pdata.2c
Pdata.2c$Var1 <- factor(Pdata.2c$Var1, levels = samp_order)
Pdata.2c$Var2 <- factor(Pdata.2c$Var2, levels = samp_order)
Pdata.2c$Var1_lab <- stringr::str_extract(Pdata.2c$Var1,".$")
Pdata.2c$Var2_lab <- stringr::str_extract(Pdata.2c$Var2,".$")

# library(ggplot2)
ggplot(data = Pdata.2c , aes(x=Var1, y=Var2, fill=value)) + 
geom_tile() + labs(x= "" , y ="" , fill = "Pairwise\nPearson\ncorrelation")+
scale_fill_gradient(limits=c(0,1), low = "black" , high = "red"  )+
theme(axis.text = element_text(size = 18, color = "black"),
	axis.title.x = element_text(size = 20, hjust =0),
	aspect.ratio = 1,
	panel.background = element_blank(),
	plot.margin=unit(c(0.1,0.1,.7,.7),"in")) +
scale_x_discrete(breaks = Pdata.2c$Var1, labels = Pdata.2c$Var1_lab ) +
scale_y_discrete(breaks = Pdata.2c$Var2, labels = Pdata.2c$Var2_lab, limits = rev) +
geom_segment(aes(x = 0.5, xend = 13.5, y = 6.5, yend = 6.5), color = "white", linetype="11", size = .8)+
geom_segment(aes(x = 0.5, xend = 13.5, y = 9.5, yend = 9.5), color = "white", linetype="11", size = .8)+
geom_segment(aes(x = 4.5, xend = 4.5, y = 0.5, yend = 13.5), color = "white", linetype="11", size = .8)+
geom_segment(aes(x = 7.5, xend = 7.5, y = 0.5, yend = 13.5), color = "white", linetype="11", size = .8)+
geom_text(inherit.aes=F,data = data.frame(x=c(2.5,6,10.5),y=rep(-0.6,3),lab=c("A3A","A3B-CTD","A3B-full")), aes(x=x,y=y,label=lab), size = 10)+
geom_text(inherit.aes=F,data = data.frame(y=c(11.5,8,3.5),x=rep(-0.7,3),lab=c("A3A","A3B-CTD","A3B-full")), aes(x=x,y=y,label=lab), angle=90, size = 10)+
coord_cartesian(xlim= c(1,13), ylim=c(1,13),clip = "off") -> g.2c

ggsave(filename= paste0(outdir,"Fig2C_sample_correlations.pdf"),plot=g.2c, height = 10, width = 10, units = "in" )

####################################################################################
### Figure 2D
## Uracilation in hairpins by looplen and looppos

### statistical analysis for Fig 2D
filter(X , U < 0.8 ) %>%  
left_join( select( filter(POSData, ss >= 12, looplen %in% 3:8 & looppos > 0 & looppos <= looplen) , POS,  looplen, looppos,ss) , .) %>%  
left_join(.,Samples, by="samp") %>%
summarySE( measurevar="U", groupvars=c("stype","samp","looplen","looppos")) %>% 
mutate(U=U*1e3) -> Data.2d
### one-waye ANOVA to determine if there is a significant difference between the samples in each looplen-looppos group
select(Data.2d,looplen, looppos) %>% distinct() -> pvals_2d 
for(i in 1:nrow(pvals_2d)){
LL <- pvals_2d$looplen[i]
LP <- pvals_2d$looppos[i]
filter(Data.2d, looplen == LL, looppos == LP ) %>%
oneway.test(data=.,U ~ stype, var.equal = FALSE) %$% p.value -> pvals_2d$pval_ANNOVA[i]
}
### multiple-hypothesis p-value adjustment and filtering for significant looplen-looppos groups
p.adjust(pvals_2d$pval_ANNOVA, method = "BH") -> pvals_2d$pval_ANNOVA_adj
filter(pvals_2d, pval_ANNOVA_adj < 0.05) %>% select(looplen , looppos) -> pvals_2d_sig
### post-hoc test for the significant looplen-looppos groups
for(i in 1:nrow(pvals_2d_sig)){
LL <- pvals_2d_sig$looplen[i]
LP <- pvals_2d_sig$looppos[i]
	filter(Data.2d, looplen == LL, looppos == LP, stype %in% c("A3A","EV")) %>%
  t.test(data=., U~stype, var.equal = FALSE) %$% p.value -> pvals_2d_sig$A3A[i]

  filter(Data.2d, looplen == LL, looppos == LP, stype %in% c("A3B_CTD","EV")) %>%
  t.test(data=., U~stype, var.equal = FALSE) %$% p.value -> pvals_2d_sig$A3B_CTD[i]

  filter(Data.2d, looplen == LL, looppos == LP, stype %in% c("A3B_full","EV")) %>%
  t.test(data=., U~stype, var.equal = FALSE) %$% p.value -> pvals_2d_sig$A3B_full[i]
}
###
gather(pvals_2d_sig, key="stype", value="pval", A3A, A3B_CTD, A3B_full) %>%
### multiple-hypothesis p-value adjustment 
mutate(p.adj = p.adjust(pval, method = "BH")) %>%
### create labels
mutate(label = pval2asterisk(p.adj)) -> pvals_2d_sig
pvals_2d_sig$label[pvals_2d_sig$label=="n.s."] <- NA
#####################
### Creating Figure 2D
filter(X , U < 0.8, samp %in% Samples$samp[Samples$stype !="EV"] ) %>%  
left_join( select( filter(POSData, ss >= 12, looplen %in% 3:8 & looppos > 0 & looppos <= looplen) , POS,  looplen, looppos,ss) , .) %>%  
left_join(.,Samples, by="samp") %>%
summarySE( measurevar="U", groupvars=c("stype","samp","looplen","looppos")) %>% 
mutate(U=U*1e3) %>% 
summarySE(  measurevar="U", groupvars=c( "stype","looplen","looppos")) -> Pdata.2d
left_join(Pdata.2d, pvals_2d_sig, by=c("looplen","looppos","stype"))  -> Pdata.2d
######
mutate(Pdata.2d,stype= factor(stype, levels=c("A3A","A3B_CTD", "A3B_full","EV")),
		rel_looppos = c((looppos - 1) / (looplen -1) )) %>% 
ggplot() +
geom_bar(aes(x=looplen, y=U, group = as.factor(looppos), fill= rel_looppos),alpha=1,lwd = 0.2 , col="black", stat= "identity", position=position_dodge())+
theme_bw()+
facet_wrap(~stype,ncol=1,strip.position="right")+
scale_fill_gradientn(name = "Position of \nC within loop", 
breaks=c(0,1), labels=c("5'","3'"),colours = c("#39c3b6","#b9d573","#ffbd14","#f04607")) +
labs(title="Hairpin size preference (ss>=12)", x="Hairpin Loop Length",y="Uracilation Index")+
scale_y_continuous(expand = c(0,0.05), breaks=seq(0,12,3))+
scale_x_continuous(expand = c(0.01,0.01), breaks=seq(0,8))+
geom_errorbar(inherit.aes = F , mapping = aes(x = looplen ,group = rel_looppos,ymin=U-sd, ymax=U+sd), width=.2, position=position_dodge(.9)) +
geom_text(  aes(x=looplen, y=U+sd+0.25, group = as.factor(looppos),label=label), position=position_dodge(.9),size = 7) +
theme( plot.title = element_text(size = 25), 
axis.text.x = element_text(color="black",size = 22), 
axis.text.y  = element_text(size = 22), 
axis.title  = element_text(size=24) , 
strip.background = element_blank (), 
strip.text = element_text(size=24),
panel.grid.minor = element_blank(),
panel.background = element_blank(),
legend.title = element_text(size = 20,hjust=0.5), 
legend.text = element_text(size = 25) , 
legend.direction = "horizontal",
legend.position =   c(0.85,0.88)) +
coord_cartesian(ylim=c(0,14.2))+
guides(fill = guide_colourbar(barwidth = 10, barheight = 1, label.position = "bottom", title.position = "top")) -> g.2d
ggsave(filename= paste0(outdir,"Fig2D_Looplen_Looppos.pdf"),plot=g.2d, height = 10, width = 10, units = "in" )

####################################################################################
### Figure 2E
## Uracilation in hairpins by loop sequence

filter(POSData , looplen %in% 3:5 & looppos >0 & looppos <= looplen & ss >=10) -> POSX
POSX %<>% mutate(loopseq = get_hairpin_seq(.,onlyLoop=T)) 

# count rows by group
group_by(POSX,looplen,looppos,loopseq) %>% 
dplyr::summarise(n_total = n()) -> SX 

select(POSX, POS, looplen, looppos, loopseq) %>% 
left_join(filter(X , U < 0.8))  %>% 
left_join(.,Samples, by="samp") -> POSX_U
##################################################
#==> Export Data for Supp. Table 3 & Table 1 
###### Supp. Table 3
filter(POSX_U, stype != "EV") %>% 
summarySE(., measurevar="U", groupvars=c("stype","samp_rep","loopseq")) %>% 
mutate( UI = U*1e3) %>% 
summarySE(., measurevar="UI", groupvars=c("stype","loopseq")) %>% 
left_join(.,SX, by="loopseq") %>% filter(n_total>=5) -> LOOPS_UI
write.table(LOOPS_UI, paste0(topdir,"Hairpin_LoopSeq_UI_ss10.ST3.txt"),row.names=F,col.names=T,quote=F,sep="\t")
###### Table 1
# top pentaloops A3A
Hairpin_LoopSeq_UI %>% select(stype,loopseq,UI,looplen,looppos,n_total ) %>% 
filter(looplen==5, stype !="EV") %>% 
pivot_wider(names_from="stype", values_from="UI") %>%
slice_max(n=20,A3A)  
# # A tibble: 20 × 7
#    loopseq   looplen looppos n_total   A3A A3B_CTD A3B_full
#    <chr>       <int>   <int>   <int> <dbl>   <dbl>    <dbl>
#  1 (TCT.C.G)       5       4       5 25.6    0.427     1.92
#  2 (TAT.C.G)       5       4      29 15.4    2.80      6.82
#  3 (TGT.C.G)       5       4      19 13.9    3.30      7.38
#  4 (GCT.C.A)       5       4      16 13.1    1.58      2.17
#  5 (TAT.C.T)       5       4      26 12.6    4.88     12.7 
#  6 (GAT.C.G)       5       4      27 12.2    2.55      4.67
#  7 (CCT.C.A)       5       4       7 11.6    1.92      1.81
#  8 (GCT.C.T)       5       4      13 10.0    1.80      2.57
#  9 (GTT.C.G)       5       4      23  8.67   1.56      3.59
# 10 (TTT.C.G)       5       4      30  7.91   1.30      5.13
# 11 (GAT.C.A)       5       4      31  7.83   0.873     3.08
# 12 (AAT.C.A)       5       4      41  5.71   1.68      3.08
# 13 (CGT.C.A)       5       4      28  5.51   2.11      4.21
# 14 (AAT.C.G)       5       4      25  5.32   2.41      3.88
# 15 (TT.C.GT)       5       3      18  5.09   2.04      5.19
# 16 (TTT.C.T)       5       4      38  4.87   1.80      4.37
# 17 (CATT.C.)       5       5      22  4.83   8.92     20.6 
# 18 (GTT.C.A)       5       4      33  4.81   0.318     2.35
# 19 (GAT.C.T)       5       4      28  4.64   1.07      2.96
# 20 (ATT.C.G)       5       4      16  4.59   0.954     3.72

# top pentaloops A3B
LOOPS_UI %>% select(stype,loopseq,UI,looplen,looppos,n_total) %>% 
filter(looplen==5) %>% 
pivot_wider(names_from="stype", values_from="UI") %>% 
slice_max(n=20,A3B_full) 
# # A tibble: 20 × 7
#    loopseq   looplen looppos n_total    A3A A3B_CTD A3B_full
#    <chr>       <int>   <int>   <int>  <dbl>   <dbl>    <dbl>
#  1 (CATT.C.)       5       5      22  4.83    8.92     20.6 
#  2 (TCAT.C.)       5       5      34  1.14    9.21     18.7 
#  3 (TTAT.C.)       5       5      46  2.66    5.41     13.5 
#  4 (CCAT.C.)       5       5      21  1.35    5.47     13.0 
#  5 (TAT.C.T)       5       4      26 12.6     4.88     12.7 
#  6 (ATTT.C.)       5       5      44  1.74    5.20     12.4 
#  7 (CAT.C.C)       5       4      20  1.70    2.54     11.5 
#  8 (CAAT.C.)       5       5      23  1.74    6.58     11.2 
#  9 (AAAT.C.)       5       5      38  0.926   4.24     11.0 
# 10 (TTGT.C.)       5       5      24  1.57    2.29     10.5 
# 11 (CGGT.C.)       5       5      18  0.809   3.47      9.72
# 12 (CTAT.C.)       5       5      15  0.926   5.41      9.70
# 13 (ACAT.C.)       5       5      20  0.924   5.13      9.62
# 14 (CCGT.C.)       5       5      19  0.419   5.52      9.15
# 15 (CGAT.C.)       5       5      27  1.00    2.94      9.07
# 16 (GGT.C.T)       5       4       9  0.484   2.37      8.47
# 17 (GT.C.CG)       5       3      13  0.840   0.372     8.30
# 18 (ATAT.C.)       5       5      37  0.821   4.98      8.23
# 19 (TAAT.C.)       5       5      23  1.34    4.87      8.19
# 20 (TGT.C.G)       5       4      19 13.9     3.30      7.38

#########################
summarySE(POSX_U, measurevar="U", groupvars=c("stype","loopseq")) %>% 
mutate(U=U*1e3) %>% 
select(1:4) %>% 
left_join(.,SX) %>% 
pivot_wider(names_from = stype, values_from = c("U","N")) -> Pdata.2ef

filter(Pdata.2ef, n_total >= 5) %>% 
arrange(sample(1:nrow(.))) %>% 
mutate(n_total = ifelse(n_total >100, 100, n_total)) %>%
mutate(rel_looppos =  c((looppos - 1) / (looplen -1) )) %>% 
ggplot(aes(x=U_A3B_CTD, y= U_A3B_full, color=rel_looppos, shape = as.factor(looplen), size= n_total )) +
geom_point( aes(alpha = log(U_A3B_CTD+U_A3B_full))) +
scale_size(breaks=c(25,50,75,100), labels= c("25","50","75","100+")) +
geom_smooth(inherit.aes=F, aes(x=U_A3B_CTD, y= U_A3B_full),method="lm",linetype ="dashed", se = T, color = "grey", size = 1)+
theme_classic()+
scale_colour_gradientn(name = "Position of \nC within loop", 
breaks=c(0,1), labels=c("5'","3'"),colours = c("#39c3b6","#b9d573","#ffbd14","#f04607")) +
guides(alpha = "none" ,
color = guide_colourbar(barwidth = 5, barheight = 1, label.position = "bottom", title.position = "top",direction = "horizontal")) +
labs(title = "Uracilation Index",x="A3B (CTD)",y="A3B (full)",size = "Number of\noccurences in\nthe genome", shape = "Loop Size") +
theme(aspect.ratio=1,
      legend.text= element_text(size=18),
      legend.title= element_text(size = 20),
      plot.title = element_text(size = 24),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 18, color = "black")) +
scale_x_continuous(expand =c(0.02,0.02))+
scale_y_continuous(expand =c(0.02,0.02)) +
scale_shape_manual(breaks=c(3,4,5), labels= c("3","4","5"), values= c(17,15,16)) +
guides(shape = guide_legend( override.aes = list(size=5)),
       size = guide_legend(shape = 17,override.aes = list(size=c(1,3,5,7))))  -> g.2e

################### get adj R2 and p-value
filter(Pdata.2ef, looplen %in% 3:5) %>% 
filter(n_total >= 5) %>%
arrange(sample(1:nrow(.))) %$% 
lm(data =., U_A3B_full ~ U_A3B_CTD) -> lm.out
summary(lm.out)$adj.r.squared %>% round(3) ->  adjr2
capture.output (summary (lm.out)) %>%
str_subset("p-value:") %>% str_extract("p-value:\ .*") -> pval 
################### plot with adj R2 
g.2e +
annotate("text", x = 1, y = 35, label =  paste0("R^2: ",adjr2),
hjust=0,  size =6, color = "black", parse = T) +
coord_cartesian(xlim=c(0,20), ylim=c(0,40)) -> g.2e_R2
ggsave(paste0(outdir,"Fig2E_Loopseq.pdf"), g.2e_R2, height = 6.5, width = 8, units = "in")
###################

####################################################################################
### Figure 2F
###################
filter(Pdata.2ef, looplen %in% 3:5) %>% 
filter(n_total >= 5) %>%
arrange(sample(1:nrow(.))) %>% 
mutate(n_total = ifelse(n_total >100, 100, n_total)) %>%
mutate(rel_looppos =  c((looppos - 1) / (looplen -1) )) %>% 
ggplot(aes(x=U_A3A, y= U_A3B_full, color=rel_looppos, shape = as.factor(looplen), size= n_total )) +
geom_point( aes(alpha = log(U_A3A+U_A3B_full))) +
scale_size(breaks=c(25,50,75,100), labels= c("25","50","75","100+")) +
theme_classic()+
scale_colour_gradientn(name = "Position of \nC within loop", 
breaks=c(0,1), labels=c("5'","3'"),colours = c("#39c3b6","#b9d573","#ffbd14","#f04607")) +
guides(alpha = "none" ,
color = guide_colourbar(barwidth = 5, barheight = 1, label.position = "bottom", title.position = "top",direction = "horizontal")) +
labs(title = "Uracilation Index",x="A3A",y="A3B",size = "Number of\noccurences in\nthe genome", shape = "Loop Size") +
theme(aspect.ratio=1,
	legend.text= element_text(size=18),
	legend.title= element_text(size = 20),
	plot.title = element_text(size = 24),
	axis.title = element_text(size = 22),
	axis.text = element_text(size = 18, color = "black")) +
scale_x_continuous(expand =c(0.02,0.02))+
scale_y_continuous(expand =c(0.02,0.02)) +
scale_shape_manual(breaks=c(3,4,5), labels= c("3","4","5"), values= c(17,15,16)) +
guides(shape = guide_legend( override.aes = list(size=5)),
	size = guide_legend(shape = 17,override.aes = list(size=c(1,3,5,7)))) -> g.2f
##################
# function for underlining the C in the loopseq annotation
underline_fun <- function(x){
  x %<>% stringr::str_remove_all("\\(|\\)")
  a <- stringr::str_remove_all(x ,"\\.C\\..*")
  a <- ifelse(a=="",a,paste0(a,"*"))
  b <- stringr::str_remove_all(x, ".*\\.C\\.")
  b <- ifelse(b=="",b,paste0("*",b))
  paste0(a,"underline(C)",b)
}
###################
################### get adj R2 and p-value
filter(Pdata.2ef, looplen %in% 3:5) %>% 
filter(n_total >= 5) %>%
arrange(sample(1:nrow(.))) %$% 
lm(data =., U_A3A ~ U_A3B_full ) -> lm.out
summary(lm.out)$adj.r.squared %>% round(3) ->  adjr2
###################
g.2f +
ggrepel::geom_text_repel(data= filter(Pdata.2ef, looplen %in% 3:5 , n_total >= 5, U_A3A >8 | U_A3B_full > 12),
aes( label = underline_fun(loopseq) ,size =(U_A3A^1.1 + U_A3B_full^1.1) ),
show.legend = F, nudge_x = 0.1, nudge_y = 0.1,colour = "#670000", parse=TRUE )+
annotate("text", x = 18, y = 35, label =  paste0("R^2: ",adjr2),
hjust=0,  size =6, color = "black", parse = T) -> g.2f_Loopseq
ggsave(paste0(outdir,"Fig2F_loopseq.pdf"), g.2f_Loopseq, height = 6.5, width = 8, units = "in")

####################################################################################
####################################################################################
### Figure 3
### River Plots
load(paste0(topdir,"Uracilation_Data.RData"))
ref_genome <- "BSgenome.Ecoli.WayneState.BH214V4"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns);library(RamPack)
library(magrittr);library(ggplot2)
###########
# function to convert X data.frame to GRanges
Xdata2GR <- function(X){
mutate(X,chr="BH214V4",start=POS,end=POS, strand = c("+","-")[match(strand,c(1,-1))],REF=ifelse(ref==2,"C","G"), ALT=ifelse(ref==2,"T","A")) %>% 
select(chr,start,end,strand,ssGroup,looplen,looppos,LDLGST,stype,REF,ALT,U) %>%
GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE,seqinfo= seqinfo(BH214)) %>%
{split(.,.$stype)} 
}
###########
## Extract the data for all positions
th <- 0.04
dplyr::filter(POSData ) %>% 
select("POS","strand","ref","ss","ssGroup","minus0","plus1","looplen","looppos","LDLGST") -> xposdata
filter(X ,U >= th & U < 0.8, samp %in% Samples$samp[Samples$stype != "EV" ] ) %>% 
left_join(xposdata ,) %>%  
left_join(.,Samples, by="samp") ->  D001

##################################################################
########### ALL POSITIONS
##### Normalization Matrix ( counting all positions)
PD_all <- list()
xposdata %>% 
# we need to remove the first and last 4 positions of the genome to avoid an error
filter(POS >= 4 & POS <= length(BH214$BH214V4)-4 ) %>% 
mutate(chr="BH214V4",start=POS,end=POS, strand = c("+","-")[match(strand,c(1,-1))],REF=ifelse(ref==2,"C","G"), ALT=ifelse(ref==2,"T","A")) %>% 
select(chr,start,end,strand,ssGroup,looplen,looppos,LDLGST,REF,ALT) %>%
GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE,seqinfo= seqinfo(BH214)) %>%
get_mut_context_matrix( . , context = "NNNCNNN" ,custom_name="genome") -> PD_all[["norm_mat"]] 
rownames(PD_all[["norm_mat"]])  %<>% stringr::str_replace_all('\\[C>T\\]',"[C>U]")
#####
filter(D001, POS >= 4 & POS <= length(BH214$BH214V4)-4 ) %>%
Xdata2GR() %>%
get_mut_context_matrix( . , context = "NNNCNNN" ) -> PD_all[["plot_data"]]
rownames(PD_all[["plot_data"]]) %<>% stringr::str_replace_all('\\[C>T\\]',"[C>U]") 
##### normalize the data
match(rownames(PD_all[["plot_data"]][,-2]),
stringr::str_replace_all(rownames(PD_all[["norm_mat"]]),'\\[C>T\\]',"[C>U]")) -> rowid
PD_all[["plot_data"]] / (PD_all[["norm_mat"]][rowid, ]) -> PD_all[["plot_data_norm"]] 
#### get Facet Names
paste0(colnames(PD_all[["plot_data_norm"]][,-2])," (n=",
colSums(PD_all[["plot_data"]][,-2]),")") -> fname
#####
RivPlot(PD_all[["plot_data_norm"]][,-2], title = "All Genomic Positions",fontsize=4,is_norm=TRUE,facet=fname) +
labs(y="Normalized Rate") -> g.3a
ggsave(filename= paste0(outdir,"Fig3A_All_newcolors_normalized.pdf"),plot=g.3a , height = 6, width = 8, units = "in" )
####################################################################################################################################
########### ALL HAIRPIN POSITIONS
PD_hp <- list()
##### Normalization Matrix ( counting all hairpin positions)
xposdata %>% 
# we need to remove the first and last 4 positions of the genome to avoid an error
filter(POS >= 4 & POS <= length(BH214$BH214V4)-4 ) %>% 
filter(ss>=10 , looplen %in% 3:6, looppos > 0 , looppos <=looplen) %>% # filter for hairpins 
mutate(chr="BH214V4",start=POS,end=POS, strand = c("+","-")[match(strand,c(1,-1))],REF=ifelse(ref==2,"C","G"), ALT=ifelse(ref==2,"T","A")) %>% 
select(chr,start,end,strand,ssGroup,looplen,looppos,LDLGST,REF,ALT) %>%
GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE,seqinfo= seqinfo(BH214)) %>%
get_mut_context_matrix( . , context = "NNNCNNN" ,custom_name="genome") -> PD_hp[["norm_mat"]] 
rownames(PD_hp[["norm_mat"]])  %<>% stringr::str_replace_all('\\[C>T\\]',"[C>U]")
#####
filter(D001, POS >= 4 & POS <= length(BH214$BH214V4)-4 ) %>%
filter(ss>=10 , looplen %in% 3:6, looppos > 0 , looppos <=looplen) %>% # filter for hairpins 
Xdata2GR() %>%
get_mut_context_matrix( . , context = "NNNCNNN" ) -> PD_hp[["plot_data"]]
rownames(PD_hp[["plot_data"]]) %<>% stringr::str_replace_all('\\[C>T\\]',"[C>U]") 
##### normalize the data
match(rownames(PD_hp[["plot_data"]][,-2]),
stringr::str_replace_all(rownames(PD_hp[["norm_mat"]]),'\\[C>T\\]',"[C>U]")) -> rowid
PD_hp[["plot_data"]] / (PD_hp[["norm_mat"]][rowid, ]) -> PD_hp[["plot_data_norm"]] 
#### get Facet Names
paste0(colnames(PD_hp[["plot_data_norm"]][,-2])," (n=",
colSums(PD_hp[["plot_data"]][,-2]),")") -> fname
#####
RivPlot(PD_hp[["plot_data_norm"]][,-2], title = "All Genomic Positions",fontsize=4,is_norm=TRUE,facet=fname) +
labs(y="Normalized Rate") -> g.3b
ggsave(filename= paste0(outdir,"Fig3B_Hairpins_newcolors_normalized.pdf"),plot=g.3b , height = 6, width = 8, units = "in" )
####################################################################################################################################
###########  (4/4) HAIRPIN POSITIONS
PD_hp44 <- list()
##### Normalization Matrix ( counting the 4/4 hairpin positions)
xposdata %>% 
# we need to remove the first and last 4 positions of the genome to avoid an error
filter(POS >= 4 & POS <= length(BH214$BH214V4)-4 ) %>% 
filter(ss>=10 , looplen %in% 3:6, looppos > 0 , looppos <=looplen) %>% # filter for hairpins 
filter(looplen==4,looppos==4) %>%
mutate(chr="BH214V4",start=POS,end=POS, strand = c("+","-")[match(strand,c(1,-1))],REF=ifelse(ref==2,"C","G"), ALT=ifelse(ref==2,"T","A")) %>% 
select(chr,start,end,strand,ssGroup,looplen,looppos,LDLGST,REF,ALT) %>%
GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE,seqinfo= seqinfo(BH214)) %>%
get_mut_context_matrix( . , context = "NNNC" ,custom_name="genome") -> PD_hp44[["norm_mat"]] 
rownames(PD_hp44[["norm_mat"]])  %<>% stringr::str_replace_all('\\[C>T\\]',"[C>U]")
#####
filter(D001, POS >= 4 & POS <= length(BH214$BH214V4)-4 ) %>%
filter(ss>=10 , looplen %in% 3:6, looppos > 0 , looppos <=looplen) %>% # filter for hairpins 
filter(looplen==4,looppos==4) %>%
Xdata2GR() %>%
get_mut_context_matrix( . , context = "NNNC" ) -> PD_hp44[["plot_data"]]
rownames(PD_hp44[["plot_data"]]) %<>% stringr::str_replace_all('\\[C>T\\]',"[C>U]") 
##### normalize the data
match(rownames(PD_hp44[["plot_data"]][,-2]),
stringr::str_replace_all(rownames(PD_hp44[["norm_mat"]]),'\\[C>T\\]',"[C>U]")) -> rowid
PD_hp44[["plot_data"]] / (PD_hp44[["norm_mat"]][rowid, ]) -> PD_hp44[["plot_data_norm"]] 
#### get Facet Names
paste0(colnames(PD_hp44[["plot_data_norm"]][,-2])," (n=",
colSums(PD_hp44[["plot_data"]][,-2]),")") -> fname
#####
RivPlot(PD_hp44[["plot_data_norm"]][,-2], title = "All Genomic Positions",fontsize=4,is_norm=TRUE,facet=fname) +
scale_x_discrete(expand=c(0.05,0.05))+
labs(y="Normalized Rate") -> g.3c
ggsave(filename= paste0(outdir,"Fig3C_HP44_newcolors_normalized.pdf"),plot=g.3c, height = 6, width = 8, units = "in" )
####################################################################################################################################
###########  (5/5) HAIRPIN POSITIONS
PD_hp55 <- list()
##### Normalization Matrix ( counting the 5/5 hairpin positions)
xposdata %>% 
# we need to remove the first and last 4 positions of the genome to avoid an error
filter(POS >= 4 & POS <= length(BH214$BH214V4)-4 ) %>% 
filter(ss>=10 , looplen %in% 3:6, looppos > 0 , looppos <=looplen) %>% # filter for hairpins 
filter(looplen==5,looppos==5) %>%
mutate(chr="BH214V4",start=POS,end=POS, strand = c("+","-")[match(strand,c(1,-1))],REF=ifelse(ref==2,"C","G"), ALT=ifelse(ref==2,"T","A")) %>% 
select(chr,start,end,strand,ssGroup,looplen,looppos,LDLGST,REF,ALT) %>%
GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE,seqinfo= seqinfo(BH214)) %>%
get_mut_context_matrix( . , context = "NNNNC" ,custom_name="genome") -> PD_hp55[["norm_mat"]] 
rownames(PD_hp55[["norm_mat"]])  %<>% stringr::str_replace_all('\\[C>T\\]',"[C>U]")
#####
filter(D001, POS >= 4 & POS <= length(BH214$BH214V4)-4 ) %>%
filter(ss>=10 , looplen %in% 3:6, looppos > 0 , looppos <=looplen) %>% # filter for hairpins 
filter(looplen==5,looppos==5) %>%
Xdata2GR() %>%
get_mut_context_matrix( . , context = "NNNNC" ) -> PD_hp55[["plot_data"]]
rownames(PD_hp55[["plot_data"]]) %<>% stringr::str_replace_all('\\[C>T\\]',"[C>U]") 
##### normalize the data
match(rownames(PD_hp55[["plot_data"]][,-2]),
stringr::str_replace_all(rownames(PD_hp55[["norm_mat"]]),'\\[C>T\\]',"[C>U]")) -> rowid
PD_hp55[["plot_data"]] / (PD_hp55[["norm_mat"]][rowid, ]) -> PD_hp55[["plot_data_norm"]] 
#### get Facet Names
paste0(colnames(PD_hp55[["plot_data_norm"]][,-2])," (n=",
colSums(PD_hp55[["plot_data"]][,-2]),")") -> fname
#####
RivPlot(PD_hp55[["plot_data_norm"]][,-2], title = "All Genomic Positions",fontsize=4,is_norm=TRUE,facet=fname) +
scale_x_discrete(expand=c(0.05,0.05))+
labs(y="Normalized Rate") -> g.3d
ggsave(filename= paste0(outdir,"Fig3D_HP55_newcolors_normalized.pdf"),plot=g.3d, height = 6, width = 8, units = "in" )

####################################################################################
####################################################################################
### Figure 4

A3A_colors <- c("#1b4499", "#82aefd") #  A3A colors
A3B_colors <- c("#f6530d", "#feb88d") #  A3B colors

### 4L HPL A3A
data.frame(Sample = factor(c("(TGT.C.)","(TTT.C.)","(GTT.C.)"),levels=c("(TGT.C.)","(TTT.C.)","(GTT.C.)")),
Activity = c(2.9,64.3,85.7)) -> data
mapacross(data$Sample,LOOPS_UI$loopseq[LOOPS_UI$stype=="A3A"],LOOPS_UI$UI[LOOPS_UI$stype=="A3A"]) -> data$Uracilation
data -> data1
normR = max(data$Activity)/max(data$Uracilation)
data %<>% mutate(Uracilation = Uracilation*normR) %>%
  gather(key = "var", value = "vals", -Sample)
plot_title <- "4L HPL"


ggplot(data,aes(x = Sample, y = vals, fill = var)) +
geom_col(position = "dodge", color = "black", show.legend = TRUE) +
scale_fill_manual(values = A3A_colors) +
scale_y_continuous(expand = c(0, 0), breaks= seq(0, max(data$vals) *0.95, by = 20),
sec.axis = sec_axis(trans = ~./normR  , name = "Uracilation")) +
scale_x_discrete(labels = function(x) parse(text = underline_fun(x)))+
labs(y = "Normalized % Activity", x = "Loop Sequence", fill = "",title = plot_title) +
theme_classic() +
coord_cartesian(ylim=c(0,max(data$vals)*1.1))+
theme(
  text = element_text(face = "bold", size = 40),
  axis.text.x = element_text(angle = 0, hjust = 0.4, vjust = 0.5),
  axis.title.y.left = element_text(size = 40, margin = margin(t = 1, r = 0.5, b = 0, l = 1, unit = "cm")),
  axis.title.x = element_text(size = 40, margin = margin(t = 0.7, r = 1, b = 0, l = 1, unit = "cm"), hjust = 0.5, vjust = 0.6),
  axis.title.y.right = element_text(size = 40, margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm"), vjust = 4), 
  legend.position = "top",
  legend.background = element_blank(),
  legend.margin = margin(t = 1, r = 0, b = 1, l = 0, unit = "cm"),  # Adjust the margin to move the legend upwards
  plot.title = element_text(hjust = 0.5, vjust = 2.5, margin = margin(t = 1, r = 5, b = 0, l = 0, unit = "cm")),
  plot.margin = margin(t = 0, r = 5, b = 0, l = 0, unit = "cm")) +
guides(fill=guide_legend(keywidth=2)) ->g.4cl

### 4L HPL A3B-CTD
data.frame(Sample = factor(c("(TGT.C.)","(TTT.C.)","(GTT.C.)"),levels=c("(TGT.C.)","(TTT.C.)","(GTT.C.)")),
Activity = c(98.47,100,29.7)) -> data
mapacross(data$Sample,LOOPS_UI$loopseq[LOOPS_UI$stype=="A3B_CTD"],LOOPS_UI$UI[LOOPS_UI$stype=="A3B_CTD"]) -> data$Uracilation
data -> data2
normR = max(data$Activity)/max(data$Uracilation)
data %<>% mutate(Uracilation = Uracilation*normR) %>%
  gather(key = "var", value = "vals", -Sample)
plot_title = "4L HPL"

ggplot(data,aes(x = Sample, y = vals, fill = var)) +
geom_col(position = "dodge", color = "black", show.legend = TRUE) +
scale_fill_manual(values = A3B_colors) +
scale_x_discrete(labels = function(x) parse(text = underline_fun(x)))+
scale_y_continuous(expand = c(0, 0), breaks= seq(0, max(data$vals)*0.9, by = 20),
sec.axis = sec_axis(trans = ~./normR  , name = "Uracilation")) +
labs(y = "Normalized % Activity", x = "Loop Sequence", fill = "",title = plot_title) +
theme_classic() +
coord_cartesian(ylim=c(0,max(data$vals)*1.1))+
theme(
  text = element_text(face = "bold", size = 40),
  axis.text.x = element_text(angle = 0, hjust = 0.4, vjust = 0.5),
  axis.title.y.left = element_text(size = 40, margin = margin(t = 1, r = 0.5, b = 0, l = 1, unit = "cm")),
  axis.title.x = element_text(size = 40, margin = margin(t = 0.7, r = 1, b = 0, l = 1, unit = "cm"), hjust = 0.5, vjust = 0.6),
  axis.title.y.right = element_text(size = 40, margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm"), vjust = 4), 
  legend.position = "top",
  legend.background = element_blank(),
  legend.margin = margin(t = 1, r = 0, b = 1, l = 0, unit = "cm"),  # Adjust the margin to move the legend upwards
  plot.title = element_text(hjust = 0.5, vjust = 2.5, margin = margin(t = 1, r = 5, b = 0, l = 0, unit = "cm")),
  plot.margin = margin(t = 0, r = 5, b = 0, l = 0, unit = "cm")) +
guides(fill=guide_legend(keywidth=2)) ->g.4cr

### 3L HPL A3B-CTD
data.frame(Sample = factor(c("(TT.C.)","(AT.C.)","(GT.C.)"),levels=c("(TT.C.)","(AT.C.)","(GT.C.)")),
Activity = c(47.3,74.8,58)) -> data
mapacross(data$Sample,LOOPS_UI$loopseq[LOOPS_UI$stype=="A3B_CTD"],LOOPS_UI$UI[LOOPS_UI$stype=="A3B_CTD"]) -> data$Uracilation
data -> data3
normR = max(data$Activity)/max(data$Uracilation)
data %<>% mutate(Uracilation = Uracilation*normR) %>%
  gather(key = "var", value = "vals", -Sample)
plot_title = "3L HPL"

ggplot(data,aes(x = Sample, y = vals, fill = var)) +
geom_col(position = "dodge", color = "black", show.legend = TRUE) +
scale_fill_manual(values = A3B_colors) +
scale_x_discrete(labels = function(x) parse(text = underline_fun(x)))+
scale_y_continuous(expand = c(0, 0), breaks= seq(0, max(data$vals)*0.9, by = 20),
sec.axis = sec_axis(trans = ~./normR  , name = "Uracilation")) +
labs(y = "Normalized % Activity", x = "Loop Sequence", fill = "",title = plot_title) +
theme_classic() +
coord_cartesian(ylim=c(0,max(data$vals)*1.1))+
theme(
  text = element_text(face = "bold", size = 40),
  axis.text.x = element_text(angle = 0, hjust = 0.4, vjust = 0.5),
  axis.title.y.left = element_text(size = 40, margin = margin(t = 1, r = 0.5, b = 0, l = 1, unit = "cm")),
  axis.title.x = element_text(size = 40, margin = margin(t = 0.7, r = 1, b = 0, l = 1, unit = "cm"), hjust = 0.5, vjust = 0.6),
  axis.title.y.right = element_text(size = 40, margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm"), vjust = 4), 
  legend.position = "top",
  legend.background = element_blank(),
  legend.margin = margin(t = 1, r = 0, b = 1, l = 0, unit = "cm"),  # Adjust the margin to move the legend upwards
  plot.title = element_text(hjust = 0.5, vjust = 2.5, margin = margin(t = 1, r = 5, b = 0, l = 0, unit = "cm")),
  plot.margin = margin(t = 0, r = 5, b = 0, l = 0, unit = "cm")) +
guides(fill=guide_legend(keywidth=2)) ->g.4dl

### 5L HPL A3B-CTD
data.frame(Sample = factor(c("(TCAT.C.)","(TTGT.C.)","(TTTT.C.)"),levels=c("(TCAT.C.)","(TTGT.C.)","(TTTT.C.)")),
Activity = c(71.7,67.17, 49.6)) -> data
mapacross(data$Sample,LOOPS_UI$loopseq[LOOPS_UI$stype=="A3B_CTD"],LOOPS_UI$UI[LOOPS_UI$stype=="A3B_CTD"]) -> data$Uracilation
data -> data4
normR = max(data$Activity)/max(data$Uracilation)
data %<>% mutate(Uracilation = Uracilation*normR) %>%
  gather(key = "var", value = "vals", -Sample)
plot_title = "5L HPL"

ggplot(data, aes(x = Sample, y = vals, fill = var)) +
geom_col(position = "dodge", color = "black", show.legend = TRUE) +
scale_fill_manual(values = A3B_colors) +
scale_x_discrete(labels = function(x) parse(text = underline_fun(x)))+
scale_y_continuous(expand = c(0, 0), breaks= seq(0, max(data$vals)*0.9, by = 20),
sec.axis = sec_axis(trans = ~./normR  , name = "Uracilation")) +
labs(y = "Normalized % Activity", x = "Loop Sequence", fill = "",title = plot_title) +
theme_classic() +
coord_cartesian(ylim=c(0,max(data$vals)*1.1))+
theme(
  text = element_text(face = "bold", size = 40),
  axis.text.x = element_text(angle = 0, hjust = 0.4, vjust = 0.5),
  axis.title.y.left = element_text(size = 40, margin = margin(t = 1, r = 0.5, b = 0, l = 1, unit = "cm")),
  axis.title.x = element_text(size = 40, margin = margin(t = 0.7, r = 1, b = 0, l = 1, unit = "cm"), hjust = 0.5, vjust = 0.6),
  axis.title.y.right = element_text(size = 40, margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm"), vjust = 4), 
  legend.position = "top",
  legend.background = element_blank(),
  legend.margin = margin(t = 1, r = 0, b = 1, l = 0, unit = "cm"),  # Adjust the margin to move the legend upwards
  plot.title = element_text(hjust = 0.5, vjust = 2.5, margin = margin(t = 1, r = 5, b = 0, l = 0, unit = "cm")),
  plot.margin = margin(t = 0, r = 5, b = 0, l = 0, unit = "cm")) +
guides(fill=guide_legend(keywidth=2)) ->g.4dr

ggpubr::ggarrange(plotlist=list(g.4cl,g.4cr,g.4dl,g.4dr),ncol=2,nrow=2) -> g.4c
ggsave(paste0(outdir,"Fig4CD.pdf"),g.4c, height = 18, width = 24, units = "in" )
######
# looking at correlation between Uracilation and in vitro Activity
rbind(data1,data2,data3,data4) %>% 
{cor.test(.$Activity,.$Uracilation,method="spearman")}
#         Spearman's rank correlation rho
# data:  .$Activity and .$Uracilation
# S = 116, p-value = 0.04575
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5944056 
####################################################################################
####################################################################################
#### Figure 5 A
## Patient data used to create the plots are provided.
PX <- read.table(paste0(topdir,'PX_data.txt'),sep="\t",header=T,comment.char="")

# Number of (3/3) hairpin sites:  210708
# Number of (4/4) hairpin sites:  204274
# Number of (4/5) hairpin sites:  146121
# Number of (5/5) hairpin sites:  102179
# Normalization Ratios:
normR45 = 146121/mean(c(146121,102179));
normR55 = 102179/mean(c(146121,102179))

damp=10
PX %<>% mutate( Ratio = ((n45/normR45)+damp)/((n55/normR55)+damp))
PX %>%
ggplot(aes(x=Ratio,y=frac_apobec*100,color = ttype_clr))+
geom_hline(yintercept = 50,linetype="longdash",color="black")+
geom_vline(xintercept = 1,linetype="longdash",color="black")+
geom_point(size=0.3)+
theme_classic()+
scale_x_continuous(expand=c(0,0),breaks=seq(0.5,4.5,.5),labels=c("","1","","2","","3","","4",""))+
scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10),labels=c("0","","","","","50","","","","","100"))+
scale_color_identity()+
coord_cartesian(clip = 'off', xlim=c(0.2,4.5),ylim=c(-0.01,100))+
labs(x="Position 4/5 ratio in pentaloops",y="Percent APOBEC")+
theme(legend.position = "none",
axis.text=element_text(size=16,color="black"),
axis.title=element_text(size=18,face="bold"),
axis.title.x=element_text(margin = margin(t = 10)),
axis.title.y=element_text(margin = margin(l = 5))) -> gtmp
ggplot2::ggplotGrob(gtmp) -> gtmp
gtmp$layout$z[gtmp$layout$name=="panel"] = 17
ggplotify::as.ggplot(gtmp) -> g.5A
ggsave(paste0(outdir,"Fig5A_HP_45_ratios.pdf"),g.5A, height = 4 , width = 5, units = "in")
################
## one sample t-test with null ratio =1
t.test(PX$Ratio,mu=1,alternative="two.sided",conf.level=0.95)
#         One Sample t-test
# data:  PX$Ratio
# t = 11.318, df = 3003, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 1
# 95 percent confidence interval:
#  1.036665 1.052032
# sample estimates:
# mean of x 
#  1.044349

## for samples with percent apobec > 50%
t.test(PX$Ratio[PX$frac_apobec>0.5],mu=1,alternative="two.sided",conf.level=0.95)
#         One Sample t-test
# data:  PX$Ratio[PX$frac_apobec > 0.5]
# t = 11.071, df = 119, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 1
# 95 percent confidence interval:
#  1.509467 1.731395
# sample estimates:
# mean of x 
#  1.620431

################
#### Figure 5 B
PX %$% case_when(mostA3AB=="A3A"~"A3A-most",
		mostA3AB=="A3B"~"A3B-most",
		frac_apobec< 0.05~"APOBEC-low\n(<5%)",
		TRUE~"rest") -> PX$mostA3AB2 
####
damp=10
normR33 = 210708/mean(c(210708,204274)); normR44 = 204274/mean(c(210708,204274))
filter(PX , mostA3AB2 != "rest") %>% 
	mutate(mostA3AB2 = factor(mostA3AB2, levels = c("A3A-most","APOBEC-low\n(<5%)","A3B-most")),
      Ratio = ((n33/normR33)+damp)/((n44/normR44)+damp)) -> PX

### t-test betweeb A3A-most and A3B-most
t.test(PX$Ratio[PX$mostA3AB2=="A3A-most"],PX$Ratio[PX$mostA3AB2=="A3B-most"],alternative="two.sided",conf.level=0.95)$p.value %>%
pval2asterisk() -> pvalA3AvA3B
#         Welch Two Sample t-test
# data:  PX$Ratio[PX$mostA3AB2 == "A3A-most"] and PX$Ratio[PX$mostA3AB2 == "A3B-most"]
# t = 13.587, df = 17.274, p-value = 1.181e-10
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  1.617446 2.211239
# sample estimates:
# mean of x mean of y 
#  2.895475  0.981133 
### t-test betweeb A3A-most and APOBEC low
t.test(PX$Ratio[PX$mostA3AB2=="A3A-most"],PX$Ratio[PX$mostA3AB2=="APOBEC-low\n(<5%)"],alternative="two.sided",conf.level=0.95)$p.value %>%
pval2asterisk() -> pvalA3AvLow
#         Welch Two Sample t-test
# data:  PX$Ratio[PX$mostA3AB2 == "A3A-most"] and PX$Ratio[PX$mostA3AB2 == "APOBEC-low\n(<5%)"]
# t = 14.45, df = 13.012, p-value = 2.16e-09
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  1.599200 2.161372
# sample estimates:
# mean of x mean of y 
#  2.895475  1.015189 
### t-test betweeb A3B-most and APOBEC low
t.test(PX$Ratio[PX$mostA3AB2=="A3B-most"],PX$Ratio[PX$mostA3AB2=="APOBEC-low\n(<5%)"],alternative="two.sided",conf.level=0.95)$p.value %>%
pval2asterisk() -> pvalA3BvLow
#         Welch Two Sample t-test
# data:  PX$Ratio[PX$mostA3AB2 == "A3B-most"] and PX$Ratio[PX$mostA3AB2 == "APOBEC-low\n(<5%)"]
# t = -0.62866, df = 11.06, p-value = 0.5423
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.15321058  0.08509797
# sample estimates:
# mean of x mean of y 
#  0.981133  1.015189 

PX %>%   ggplot(aes(x=mostA3AB2,y=Ratio*rnorm(n=nrow(PX), mean=1,sd=0.02),fill=mostA3AB2))+
geom_violin(alpha=0.5, position=position_dodge(.9))+
geom_jitter(width=0.1,size=0.4,alpha=0.5)+
scale_fill_manual(values =c("#fc3a46","#9e9e9e","#5a96ff"))+
theme_classic()+
theme(legend.position = "none",
plot.title=element_text(size=18,hjust=0.5),
axis.text=element_text(size=14, colour = "black"),
axis.title=element_text(size=18, colour = "black" ))+
labs(title="",y="log2 ratio",x="",fill="") +
geom_abline(intercept=0, slope=0, linetype="dashed", color = "red") +
scale_y_continuous(limits=c(0.4,5.2),trans = "log2", breaks=c(0.5,1,2,4), labels=log2(c(0.5,1,2,4))) +
annotation_logticks(sides = "l") +
annotate("text", x = 1.5, y = 4.2, label = pvalA3AvLow, size = 7) + 
geom_segment(aes(x =.9, xend = 2.1 , y = 4.1, yend = 4.1), color = "black", linewidth = 1)+
annotate("text", x = 2.5, y = 2.5, label = pvalA3BvLow, size = 5) +
geom_segment(aes(x =1.9, xend = 3.1 , y = 2.2, yend = 2.2), color = "black", linewidth = 1) +
annotate("text", x = 2, y = 5.1, label = pvalA3AvA3B, size = 7) +
geom_segment(aes(x =.9, xend = 3.1 , y = 5, yend = 5), color = "black", linewidth = 1) -> g.5B
ggsave(paste0(outdir,"Fig5B_HP_33_44_ratios.pdf"),g.5B, height = 4.2 , width = 5, units = "in")
################


#### Figures 5 C,D
################### Data for Top Panel (5C) #############################
filter(X , U < 0.8 , samp %in% Samples$samp[Samples$stype %in% c("A3A","A3B_full")] ) %>%  
left_join( select( filter(POSData, minus0==4 ) , POS,  looplen, looppos,ss,ssGroup) , .) %>%  
left_join(.,Samples, by= "samp") -> POSXDTOP_U
#####
filter(POSXDTOP_U, ssGroup=="0-7") %>% 
drop_na() %>%
summarySE( measurevar="U", groupvars=c("stype")) %>% 
mutate(U=U*1e3) -> baseRates
##########
POSXDTOP_U %>% drop_na() %>%
summarySE( measurevar="U", groupvars=c("stype","samp","ssGroup")) %>% 
mutate(U=U*1e3) %>% 
# normalize to the ss-Group 0-7
mutate(UN = ifelse(stype=="A3A",U/baseRates$U[baseRates$stype=="A3A"], U/baseRates$U[baseRates$stype=="A3B_full"])) %>% 
summarySE(  measurevar="UN", groupvars=c( "stype","ssGroup")) -> PD.5C_UPDSEQ
######
PD.5C_UPDSEQ$stype[ PD.5C_UPDSEQ$stype=="A3B_full"] <- "A3B"
PD.5C_UPDSEQ %<>% mutate(UN = ifelse(stype=="A3A",UN * -1, UN ))

#################### Data for Bottom Panel (5D) #############################
filter(X , U < 0.8 , samp %in% Samples$samp[Samples$stype %in% c("A3A","A3B_full")] ) %>%  
left_join( select( filter(POSData, minus0==4, looplen %in% 3:6 & looppos > 0 & looppos <= looplen) , POS,  looplen, looppos,ss,ssGroup) , .) %>%  
left_join(.,Samples, by= "samp")  %>% 
drop_na() %>% 
summarySE( measurevar="U", groupvars=c("stype","samp","looplen","looppos","ssGroup")) %>% 
mutate(U=U*1e3) %>% 
mutate(UN = ifelse(stype=="A3A",U/baseRates$U[baseRates$stype=="A3A"], U/baseRates$U[baseRates$stype=="A3B_full"])) %>% # normalize to the ss-Group 0-7
summarySE(  measurevar="UN", groupvars=c( "stype","looplen","looppos", "ssGroup")) -> PD.5D_UPDSEQ
PD.5D_UPDSEQ$stype[PD.5D_UPDSEQ$stype=="A3B_full"] <- "A3B"
PD.5D_UPDSEQ %<>% mutate(UN = ifelse(stype=="A3A",UN * -1, UN ))

########################################
### Plotting
library(ggnewscale)
colB <- colorRampPalette(c("#e1e6ff","#0018ca"))(8)
colA <- colorRampPalette(c("#fee0e0","#da000e"))(8)
#####
ggplot(PD.5C_UPDSEQ[PD.5C_UPDSEQ$stype=="A3B",], aes(x=ssGroup, y=UN, group = ssGroup, fill= ssGroup))+
geom_bar( ,alpha=1,lwd = 0.2 , col="black", stat= "identity", position=position_dodge())+
scale_fill_manual(values=colB) + 
new_scale_fill() +
geom_bar( data=PD.5C_UPDSEQ[PD.5C_UPDSEQ$stype=="A3A",],aes(x=ssGroup, y=UN, group = ssGroup, fill= ssGroup),alpha=1,lwd = 0.2 , col="black", stat= "identity", position=position_dodge())+
scale_fill_manual(values=colA) +
theme_bw()+
geom_errorbar( data=PD.5C_UPDSEQ, inherit.aes = T , mapping = aes( group = ssGroup,ymin=UN-ci, ymax=UN+ci), width=.2, position=position_dodge(.9)) + 
geom_segment(aes(x=0, xend=8.8, y=0, yend=0), linetype="dashed",color="black", linewidth=1) +
annotate(geom="text", x=1, y=-1.8,size=9, label="A3A",color="red",hjust=1) +
annotate(geom="text", x=1, y=1.85,size=9, label="A3B",color="Blue",hjust=0) +
labs(title="All TpC sites", x="Stem Strength",y="" , fill= "Hairpin\nStem \nStrength")+
theme( plot.title = element_text(size = 25), 
axis.text.x = element_text(color="black",size = 22), 
axis.text.y  = element_text(size = 22), 
axis.title  = element_text(size=24) , 
strip.background = element_rect(fill="white"), 
strip.text = element_text(size=24),
legend.title = element_text(size = 20,hjust=0.5), 
legend.text = element_text(size = 22) , 
legend.direction = "horizontal",
legend.position =   "none",
plot.margin = margin(.5, 0, 0, 0, "in")) +
scale_y_continuous(expand = c(0,0), breaks=seq(-8,8,2),labels=abs(seq(-8,8,2)))+
guides(fill = guide_colourbar(barwidth = 10, barheight = 1, label.position = "bottom", title.position = "top")) +
coord_flip(ylim=c(-8,8),xlim=c(0.8,8.2),clip="off") +
annotate(geom="text", x=12, y=0,size=12, label="italic('E. coli')~ ' UPD-Seq'",parse=T,color="Black",hjust=0.5) -> g.5C_UPDSEQ
#####
ggplot(PD.5D_UPDSEQ) +
geom_bar(data = filter(PD.5D_UPDSEQ, stype=="A3B"), aes(x=looppos, y=UN, group = as.factor(ssGroup), fill= ssGroup),alpha=1,lwd = 0.2 , col="black", stat= "identity", position=position_dodge( preserve = "single"))+
scale_fill_manual(values=colB) + 
new_scale_fill() +
geom_bar(data = filter(PD.5D_UPDSEQ, stype=="A3A"), aes(x=looppos, y=UN, group = as.factor(ssGroup), fill= ssGroup),alpha=1,lwd = 0.2 , col="black", stat= "identity", position=position_dodge( preserve = "single"))+
scale_fill_manual(values=colA) +
theme_bw()+
coord_flip(ylim=c(-80,80)) +
facet_wrap(~looplen ,labeller = as_labeller(c("3"="3nt loop", "4"="4nt loop","5"="5nt loop","6"="6nt loop")), scales= "free_y", ncol=1,strip.position = "right") +
geom_segment(data= data.frame(x=.3, xend=3.7, y=0, yend=0,looplen=3),
aes(x=x, xend=xend, y=y, yend=yend), linetype="solid",color="black", size=.8,inherit.aes=FALSE) +
geom_segment(data= data.frame(x=.3, xend=4.7, y=0, yend=0,looplen=4),
aes(x=x, xend=xend, y=y, yend=yend), linetype="solid",color="black", size=.8,inherit.aes=FALSE) +
geom_segment(data= data.frame(x=.3, xend=5.7, y=0, yend=0,looplen=5),
aes(x=x, xend=xend, y=y, yend=yend), linetype="solid",color="black", size=.8,inherit.aes=FALSE) +
geom_segment(data= data.frame(x=.3, xend=6.7, y=0, yend=0,looplen=6),
aes(x=x, xend=xend, y=y, yend=yend), linetype="solid",color="black", size=.8,inherit.aes=FALSE) +
labs(title="Hairpin loops",fill="Stem Strength", x="Position of TpC within the loop",y="Normalized Uracilation Index") +
scale_y_continuous(expand = c(0,0), breaks=seq(-80,80,20),labels=abs(seq(-80,80,20)))+
scale_x_reverse(expand = c(0,0) , breaks=seq(1,6,1)) +
geom_errorbar(inherit.aes = F , mapping = aes(x = looppos ,group = ssGroup,ymin=UN-ci, ymax=UN+ci), width=.2, position=position_dodge(.9)) +
theme( plot.title = element_text(size = 25), 
axis.text.x = element_text(color="black",size = 22), 
axis.text.y  = element_text(size = 22), 
axis.title.x = element_text(color="black",size = 22), 
axis.title.y  = element_text(size=24) , 
strip.background = element_rect(fill="white"), 
strip.text = element_text(size=24),
legend.title = element_text(size = 15,hjust=0.5), 
legend.text = element_text(size = 14) , 
legend.position =  "none",
legend.background = element_rect(fill=NA))+
scale_fill_manual(values =  colA ) -> g.5D_UPDSEQ 

ggpubr::ggarrange(g.5C_UPDSEQ+theme(plot.margin = margin(0, 1, 0, 0, "in")), g.5D_UPDSEQ, ncol = 1, heights = c(2,5)) +
theme(plot.margin = unit(c(1, .2, .2, .2), "in")) -> G.5CD_UPDSEQ

########################################################################## 
##################### import mutation Data from Human tumors
#################### Data for Top Panel (5C) #############################
TumorData5C <- read.table(paste0(topdir,"Tumors_A3A_A3B_TpC_ssbins.txt"),header=T)
TumorData5C <- mutate(TumorData5C,CI=sd_TpC * 1.96 ,
					Rate = {ifelse(A3AB=="A3A", (rate_TpC * -1 ),rate_TpC )}) %>%
		mutate(range = factor(range, levels=c('0-4','5-7','8-11','12-15','16+')))
#################### Data for Bottom Panel (5D) #############################
TumorData5D <- read.table(paste0(topdir,"Tumors_A3A_A3B_TpC_ssbins_looplen_looppos.txt"),header=T)
TumorData5D <- mutate(TumorData5D, CI=sd*1.96,
		Rate = {ifelse(Sample=="A3A", (rate  * -1 ),rate )},
		range = factor(ss_range, levels=c('0-4','5-7','8-11','12-15','16+')))
### set a minimum number of mutations per bin to be included in the plot
nmin=2;
TumorData5D %<>% mutate(Rate = ifelse(nmut < nmin, NA, Rate),
						CI= ifelse(nmut < nmin , NA, CI))  
##########################
### Plotting
library(ggplot2);library(ggnewscale)
colB <- colorRampPalette(c("#e1e6ff","#0018ca"))(6)
colA <- colorRampPalette(c("#fee0e0","#da000e"))(6)
###
ggplot(TumorData5C[TumorData5C$A3AB=="A3B",], aes(x=range, y=Rate, group = range, fill= range))+
geom_bar( ,alpha=1,lwd = 0.2 , col="black", stat= "identity", position=position_dodge())+
scale_fill_manual(values=colB) + 
new_scale_fill() +
geom_bar( data=TumorData5C[TumorData5C$A3AB=="A3A",],aes(x=range, y=Rate, group = range, fill= range),alpha=1,lwd = 0.2 , col="black", stat= "identity", position=position_dodge())+
scale_fill_manual(values=colA) +
theme_bw()+
geom_errorbar( data=TumorData5C, inherit.aes = T , mapping = aes( group = range,ymin=Rate-CI, ymax=Rate+CI), width=.2, position=position_dodge(.9)) + 
geom_segment(aes(x=0, xend=5.3, y=0, yend=0), linetype="dashed",color="black", size=1) +
annotate(geom="text", x=1, y=-1.2,size=9, label="A3A-most",color="red",hjust=1) +
annotate(geom="text", x=1, y=1.2,size=9, label="A3B-most",color="Blue",hjust=0) +
labs(title="All TpC sites", x="Stem Strength",y="" , fill= "Hairpin\nStem \nStrength")+
theme( plot.title = element_text(size = 25), 
axis.text.x = element_text(color="black",size = 22), 
axis.text.y  = element_text(size = 22), 
axis.title  = element_text(size=24) , 
strip.background = element_rect(fill="white"), 
strip.text = element_text(size=24),
legend.title = element_text(size = 20,hjust=0.5), 
legend.text = element_text(size = 22) , 
legend.direction = "horizontal",
plot.margin = margin(.5, 0, 0, 0, "in"),
legend.position =   "none") +
scale_y_continuous(expand = c(0,0), breaks=seq(-3,3,1),labels=abs(seq(-3,3,1)))+
guides(fill = guide_colourbar(barwidth = 10, barheight = 1, label.position = "bottom", title.position = "top")) +
coord_flip(ylim=c(-3,3),xlim=c(0.8,5.2),clip="off") +
annotate(geom="text", x=8, y=0,size=12, label="Patient tumors ",color="Black",hjust=0.5) -> g.5C_Tumors
###################
ggplot(TumorData5D) +
geom_bar(data = filter(TumorData5D, Sample=="A3B"), aes(x=looppos, y=Rate, group = as.factor(range), fill= range),alpha=1,lwd = 0.2 , col="black", stat= "identity", position=position_dodge( preserve = "single"))+
scale_fill_manual(values=colB) + 
new_scale_fill() +
geom_bar(data = filter(TumorData5D, Sample=="A3A"), aes(x=looppos, y=Rate, group = as.factor(range), fill= range),alpha=1,lwd = 0.2 , col="black", stat= "identity", position=position_dodge( preserve = "single"))+
scale_fill_manual(values=colA) +
theme_bw()+
coord_flip(ylim=c(-22,22)) +
facet_wrap(~looplen ,labeller = as_labeller(c("3"="3nt loop", "4"="4nt loop","5"="5nt loop","6"="6nt loop")), scales= "free_y", ncol=1,strip.position = "right") +
geom_segment(data= data.frame(x=.3, xend=3.7, y=0, yend=0,looplen=3),
aes(x=x, xend=xend, y=y, yend=yend), linetype="solid",color="black", size=.8,inherit.aes=FALSE) +
geom_segment(data= data.frame(x=.3, xend=4.7, y=0, yend=0,looplen=4),
aes(x=x, xend=xend, y=y, yend=yend), linetype="solid",color="black", size=.8,inherit.aes=FALSE) +
geom_segment(data= data.frame(x=.3, xend=5.7, y=0, yend=0,looplen=5),
aes(x=x, xend=xend, y=y, yend=yend), linetype="solid",color="black", size=.8,inherit.aes=FALSE) +
geom_segment(data= data.frame(x=.3, xend=6.7, y=0, yend=0,looplen=6),
aes(x=x, xend=xend, y=y, yend=yend), linetype="solid",color="black", size=.8,inherit.aes=FALSE) +
labs(title="Hairpin loops",fill="Stem Strength", x="Position of TpC within the loop",y="Relative mutation Frequency") +
scale_y_continuous(expand = c(0,0), breaks=seq(-20,20,10),labels=c(seq(20,0,-10),10,20))+
scale_x_reverse(expand = c(0,0) , breaks=seq(1,6,1)) +
geom_errorbar(inherit.aes = F , mapping = aes(x = looppos ,group = range,ymin=Rate-CI, ymax=Rate+CI), width=.2, position=position_dodge(.9)) +
theme( plot.title = element_text(size = 25), 
axis.text.x = element_text(color="black",size = 22), 
axis.text.y  = element_text(size = 22), 
axis.title.x = element_text(color="black",size = 22), 
axis.title.y  = element_text(size=24) , 
strip.background = element_rect(fill="white"), 
strip.text = element_text(size=24),
legend.title = element_text(size = 15,hjust=0.5), 
legend.text = element_text(size = 14) , 
legend.position =  "none",
legend.background = element_rect(fill=NA))+
scale_fill_manual(values =  colA ) -> g.5D_Tumors
###################
ggpubr::ggarrange(g.5C_Tumors+theme(plot.margin = margin(0, 1, 0, 0, "in")), g.5D_Tumors, ncol = 1, heights = c(2,5)) +
  theme(plot.margin = unit(c(1, .2, .2, .2), "in")) -> G.5CD_Tumors
######################
ggpubr::ggarrange(G.5CD_UPDSEQ,G.5CD_Tumors, ncol = 2 ) +
  theme(plot.margin = unit(c(.2, .2, .2, .2), "in")) -> G.5CD
ggsave(paste0(outdir,"Fig5_CD_UPDSEQ_PatientTumors.pdf") ,  G.5CD,  height = 15 , width = 18, units = "in")


#################################
###### for supplementary figure 6
### Generating a made up list of mutations to create supplementary figure explaining the river plot
matrix(c(2,10,10,5,10,3,7,2,4,5,2,3,4), nrow = 13, ncol = 1,
 dimnames = list(c("AT[C>U]A", "TT[C>U]C", "TT[C>U]G", "TT[C>U]T", "CA[C>U]A", "CA[C>U]C", "CA[C>U]G", "GT[C>U]T", "CT[C>U]A", "CA[C>U]C", "CG[C>U]G", "CC[C>U]T","TA[C>U]A"), "Sample")) -> RivPlot_tmp

## counting the mutations in a 4x4 matrix
matrix(0,nrow=4, ncol = 4, dimnames =list(c("A","C","G","T"),1:4)) -> Seqlogo_tmp
for(i in 1:nrow(RivPlot_tmp)){
	for(j in c(1,2,4)){
		x <- c(1,2,0,8)[j]
  Seqlogo_tmp[substr(rownames(RivPlot_tmp)[i],x,x),j] <- Seqlogo_tmp[substr(rownames(RivPlot_tmp)[i],x,x),j] + RivPlot_tmp[i,1]
}}
Seqlogo_tmp["C",3] <- colSums(Seqlogo_tmp)[1]

### creaing Seqlogo and River plot
ggseqlogo( Seqlogo_tmp, method = 'bits' )+
 scale_x_continuous(breaks=1:4, labels=c("-2","-1","0","+1"))+
 labs(title = "SeqLogo", x = "Position", y = "Bits")+
   theme(axis.text = element_text(size= 18),
   plot.title = element_text(size= 22),
     axis.title = element_text(size= 20))

RivPlot(RivPlot_tmp, title = "River Plot",fontsize=8,is_norm=TRUE ) -> g.6S
ggsave(paste0(outdir,"Fig6S_RiverPlot.pdf"),g.6S, height = 4 , width = 5, units = "in") 
#################################