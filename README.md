
## Distinguishing preferences of human APOBEC3A and APOBEC3B for cytosines in hairpin loops, and reflection of these preferences in APOBEC-signature cancer genome mutations
#### Yasha Butt<sup>1,\*</sup>, Ramin Sakhtemani<sup>1,2,3,*</sup>, Rukshana Mohamad-Ramshan<sup>1</sup>, Michael S. Lawrence<sup>2,3</sup> and Ashok S. Bhagwat<sup>1,4,§</sup>

_1-_ Department of Chemistry, Wayne State University, Detroit, MI, 48202, USA<br>
_2-_ Massachusetts General Hospital Cancer Center, Boston, MA, USA<br>
_3-_ Broad Institute of MIT and Harvard, Cambridge, MA, USA<br>
_4-_ Department of Biochemistry, Microbiology and Immunology, Wayne State University School of Medicine, Detroit, MI 48201, USA<br><br>

*These authors contributed equally to this work.<br>
§ To whom correspondence should be addressed: Ashok S. Bhagwat,<br> 443 Chemistry Building,Wayne state University, Detroit, MI 48202; <br>Email- axb@chem.wayne.edu; Tel. 01-734-425-1749; Fax: (313) 577-8822.
#### [(BioRxiv: DOI: 10.1101/2023.08.01.551518)](10.1101/2023.08.01.551518)
________

### Analysis of Uracils Created by human APOBEC3A and APOBEC3B*

#### Samples.xlsx :  lists of all of the samples and their accession numbers.

**Four APOBEC3A samples are already published <br>**
A3A_F6&emsp;SRR6924522&emsp;PRJNA448166<br>
A3A_H1&emsp;SRR9864913&emsp;PRJNA448166<br>
A3A_A&emsp;SRR17822878&emsp;PRJNA801888<br>
A3A_G&emsp;SRR17822877&emsp;PRJNA801888<br>
**New samples**<br>
3 A3B-CTD + 3 Empty vector control samples<br>
6 A3B-full + 6 Empty vector control samples
________
**BH214 strain of E. _coli_ is used in the experiments.<br>The reference genome used for the analysis is BH214V4.fa available inside the BH214_genome_files directory.<br><br>
Step1: Download the raw sequences using the provided accession numbers.<br>
Step2: Sequence alignment on and extracting the depth of coverage and the nucleotide counts tables.<br>**

```
# <ref.fa> is the reference genome (BH214V.fa)
# <plasmid_seq.fa> is the plasmid sequence (A3B_plasmids.fa)
# change the <names> to match the sample names

# index fasta files
bwa index -p PAB <plasmid_seq.fa>
bwa index -p BH <ref.fa>

# alignments and filtering out plasmid reads
bwa mem -t 4 PAB <sample.fastq1> <sample.fastq2> | samtools view -b -f4 | samtools sort -n -l 6 -@ 4 - -o <sample.plasmid_aligned.bam>
bedtools bamtofastq -i <sample.plasmid_aligned.bam> -fq  <sample.NoPlasmid_R1.fq>  -fq2  <sample.NoPlasmid_R2.fq>
bwa mem -t 4 BH <sample.NoPlasmid_R1.fq> <sample..NoPlasmid_R2.fq> | samtools view -b | samtools sort -l 6 -@ 4 - -o <sample.BH214.bam>

# extracting nucleotide readcounts
bam-readcount -w0 -f <ref.fa> <sample.BH214.bam> | \
awk -F ":|\t|=" 'BEGIN {OFS = "\t"}; {print $1, $2, $3 , $4, $21 , $35, $49 , $63}' > <sample_readcount_out.txt>

# for NDC2 analysis, the depth of coverage is needed:
samtools depth -aa -m 100000 <sample.BH214.bam> ><sample.BH214.depth>
```
<br>

**Step3: Follow the instructions in Analysis.R to reproduce the analysis and create the Plots.**
________


