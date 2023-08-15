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
## then run the following code to extract the intersection of 


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
ggsave(filename= paste0(outdir,"Fig2B_HairpinSSGroup.pdf"),plot=g.2b, height = 10, width = 12, units = "in" )

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

####### Fig 2C
filter(X , U < 0.8, samp %in% Samples$samp[Samples$stype !="EV"] ) %>%  
left_join( select( filter(POSData, ss >= 12, looplen %in% 3:8 & looppos > 0 & looppos <= looplen) , POS,  looplen, looppos,ss) , .) %>%  
left_join(.,Samples, by="samp") %>%
summarySE( measurevar="U", groupvars=c("stype","samp","looplen","looppos")) %>% 
mutate(U=U*1e3) %>% 
summarySE(  measurevar="U", groupvars=c( "stype","looplen","looppos")) -> Pdata.2d
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
coord_cartesian(ylim=c(0,14))+
guides(fill = guide_colourbar(barwidth = 10, barheight = 1, label.position = "bottom", title.position = "top")) -> g.2d
ggsave(filename= paste0(outdir,"Fig2D_Looplen_Looppos.pdf"),plot=g.2d, height = 10, width = 10, units = "in" )

####################################################################################

### Figure 2E
## Uracilation in hairpins by loop sequence

filter(POSData , looplen %in% 3:6 & looppos >0 & looppos <= looplen & ss >=10) -> POSX
POSX %<>% mutate(loopseq = get_hairpin_seq(.,onlyLoop=T)) 

# count rows by group
group_by(POSX,looplen,looppos,loopseq) %>% 
dplyr::summarise(n_total = n()) -> SX 

select(POSX, POS, looplen, looppos, loopseq) %>% 
left_join(filter(X , U < 0.8))  %>% 
left_join(.,Samples, by="samp") -> POSX_U

summarySE(POSX_U, measurevar="U", groupvars=c("stype","loopseq")) %>% 
mutate(U=U*1e3) %>% 
select(1:4) %>% 
left_join(.,SX) %>% 
pivot_wider(names_from = stype, values_from = c("U","N")) -> Pdata.2ef

filter(Pdata.2ef, looplen %in% 3:5) %>% 
filter(n_total >= 5) %>%
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
g.2f +
ggrepel::geom_text_repel(data= filter(Pdata.2ef, looplen %in% 3:5 , n_total >= 5, U_A3A >8 | U_A3B_full > 12),
aes( label = underline_fun(loopseq) ,size =(U_A3A^1.1 + U_A3B_full^1.1) ),
show.legend = F, nudge_x = 0.1, nudge_y = 0.1,colour = "#670000", parse=TRUE ) -> g.2f_Loopseq
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
####################################################################################################################################
