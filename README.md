
## Distinguishing preferences of human APOBEC3A and APOBEC3B for cytosines in hairpin loops, and reflection of these preferences in APOBEC-signature cancer genome mutations
#### Yasha Butt<sup>1,\*</sup>, Ramin Sakhtemani<sup>1,2,3,*</sup>, Rukshana Mohamad-Ramshan<sup>1</sup>, Michael S. Lawrence<sup>2,3</sup> and Ashok S. Bhagwat<sup>1,4,§</sup>

_1-_ Department of Chemistry, Wayne State University, Detroit, MI, 48202, USA<br>
_2-_ Massachusetts General Hospital Cancer Center, Boston, MA, USA<br>
_3-_ Broad Institute of MIT and Harvard, Cambridge, MA, USA<br>
_4-_ Department of Biochemistry, Microbiology and Immunology, Wayne State University School of Medicine, Detroit, MI 48201, USA<br><br>

*These authors contributed equally to this work.<br>
§ To whom correspondence should be addressed: Ashok S. Bhagwat,<br> 443 Chemistry Building,Wayne state University, Detroit, MI 48202; <br>Email- axb@chem.wayne.edu; Tel. 01-734-425-1749; Fax: (313) 577-8822.<br><br>
**[(BioRxiv: DOI: 10.1101/2023.08.01.551518)](10.1101/2023.08.01.551518)**
________

## Analysis of Uracils Created by human APOBEC3A and APOBEC3B*

**Step1: Download the raw sequences using the provided accession numbers.<br>**
<h4> Samples.xlsx lists all of the samples and their accession numbers.</h4>

<h4>Four APOBEC3A samples are already published</h4>
| A3A_F6 | SRR6924522 | PRJNA448166 | <br>
| A3A_H1 | SRR9864913 | PRJNA448166 | <br>
| A3A_A  | SRR17822878 | PRJNA801888 | <br>
| A3A_G  | SRR17822877 | PRJNA801888 | <br><br>

**New samples**<br>
| 3 A3B-CTD  + 3 Empty vector control samples  | PRJNA1005650  |<br>
| 6 A3B-full (U+0009) + 6 Empty vector control samples  | PRJNA1005650  |
______
**Step2: Sequence alignment on and extracting the depth of coverage and the nucleotide counts tables.<br>**
**The reference sequences are available inside the BH214_genome_files directory.<br>**
**Follow these steps on all samples.<br><br>**


<ref.fa> is the reference genome (BH214V.fa)<br>
<plasmid_seq.fa> is the plasmid sequence (A3B_plasmids.fa)<br>
change the <names> to match the sample names<br>

\#index fasta files
```
bwa index -p PAB <plasmid_seq.fa>
bwa index -p BH <ref.fa>
```
\# alignments and filtering out plasmid reads
```
bwa mem -t 4 PAB <sample.fastq1> <sample.fastq2> | \
samtools view -b -f4 | \
samtools sort -n -l 6 -@ 4 - -o <sample.plasmid_aligned.bam>

bedtools bamtofastq -i <sample.plasmid_aligned.bam> -fq  <sample.NoPlasmid_R1.fq>  -fq2  <sample.NoPlasmid_R2.fq>
bwa mem -t 4 BH <sample.NoPlasmid_R1.fq> <sample..NoPlasmid_R2.fq> | \
samtools view -b | \
samtools sort -l 6 -@ 4 - -o <sample.BH214.bam>
```
\# extracting nucleotide readcounts
```
bam-readcount -w0 -f <ref.fa> <sample.BH214.bam> | \
awk -F ":|\t|=" 'BEGIN {OFS = "\t"}; {print $1, $2, $3 , $4, $21 , $35, $49 , $63}' > <sample_readcount_out.txt>
```
\# for NDC2 analysis, the depth of coverage is needed:
```
samtools depth -aa -m 100000 <sample.BH214.bam> ><sample.BH214.depth>
```
<br>

**Step3: Follow the instructions in Analysis.R to reproduce the analysis and create the Plots.**
________ 
<br>

## Hairpin Signature Analysis
An example MATLAB code (Hairpin_analysis.m), a test mutation set (Test_Mutations.txt), and the dependency function are available inside the Hairpin_Signature_Analysis/ directory.
____

**Hairpin Signature Analysis** is an nmf-based model to estimate the mutational activity of APOBEC3A and APOBEC3B from mutation patterns at hairpin-forming sequences.<br>
To Run the analysis on a sample's list of mutations, we first need to find the hairpin information:<br>
1- length of the hairpin loops, 2- positions of the cytosine within the loop, and 3- stem strengths are needed. <br><br>
We have created a MATLAB function to extract this information from the reference genome:<br>
_________
**get_hairpin_info(M,ref_fasta)<br><br>**
This function requires the reference fasta file, and a table of mutations with the following fields:<br>
**_ref_fasta:_ '/path/to/refernece/sequence/such/as/Homo_sapiens_assembly19.fasta'** <br>
**_M:_ A struct (table) of mutations with the following fields**<br>
_chr:_ chromosome<br>
_pos:_ position of the mutation<br> 
_ref:_ reference sequence : A/C/G/T or 1/2/3/4<br>
_alt:_ alternative sequence : A/C/G/T or 1/2/3/4
________

**Note: samtools must be available at the system level to be able to run the get_hairpin_info() function!<br><br>** 

This function Runs on a loop and, for each mutation, scans the surrounding sequence to find the potential hairpin structure.<br>
Alternatively, one can extract the needed hairpin information by scanning the entire genome for potential hairpins.<br>
**Note:** <br>
For a large number of mutations (more than ~10,000 mutations), it is likely faster to scan the entire genome for hairpins,
then map the mutations to the list of genomic positions with hairpin information and extract the loop length, loop position, and stem strength values for a set of mutations.<br>
See [https://github.com/alangenb/ApoHP](https://github.com/alangenb/ApoHP) for ApoHP: APOBEC hairpins analysis tool<br>
_____________

**Run the Hairpin Signature Analysis in MATLAB<br>**
Once the hairpin information is gathered, we can run the hairpin_signature_analysis() function.<br>
This function calculates two coefficients, hs1 and hs2, for the Hairpin Signatures HS1 (hairpin mutation pattern of APOBEC3A) and HS2 (hairpin mutation pattern of APOBEC3B).<br>

________
**hairpin_signature_analysis(Hairpins,only_tpc)<br><br>**

This function takes in (Hairpins) the list of mutations containing the hairpin information and returns the result of the hairpin signature analysis.
If a "patid"/"sampid" field is present, then the result will be repeated on each patient/sample separately. otherwise, all the mutations are assumed to be from a single sample.<br><br>
**input**<br>
**_Hairpins:_ A struct (table) of mutation/hairpin info with the following fields:** <br>  
_ref:_ reference sequence, numeric (1/2/3/4) : A=1, C=2, G=3, T=4<br>
_alt:_ alternative sequence, numeric (1/2/3/4) : A=1, C=2, G=3, T=4<br>
_looplen:_ length of the potential hairpin loop<br>
_looppos:_ position of the cytosine in the potential hairpin loop<br>
_ss:_ stem strength of the potential hairpin loop<br>
_minus0:_ base immediately 5' to the cytosine, numeric (1/2/3/4) : A=1, C=2, G=3, T=4<br>
_patid:_  (optional) patient/sample index, must be numeric<br><br>
**_only_TpC:_** true by default; if false, then it considers all of the C:G mutations, not just TpC.<br><br>
**output**<br>
A nested struct which contains:<br>
- _pat:_<br>
   - pat_id<br>
   - nC2GT: # of C>T | C>G mutations<br>
   - nTC2GT: # of C>T | C>G mutations in TpC context<br>
   - tc_frac: fraction of C:G mutations in TpC context<br>
   - hs1: coefficient of HS1 signature (APOBEC3A)<br>
   - hs2: coefficient of HS2 signature (APOBEC3B)<br>
   - log2R: log2 ratio of hs1 and hs2<br>
   - judgment: 'A3A-like' or 'A3B-like'<br>
- _counts:_ The count matrix (nx90) of mutations used in the directed NMF<br>
- _HP:_ information about the 90 hairpin groups used in the analysis and HS1 and HS2 values.<br>
________________
#### If you have any questions regarding the article or the code behind it, please contact the corresponding author or Ramin Sakhtemani rsakhtemani (at) mgh.harvard.edu
