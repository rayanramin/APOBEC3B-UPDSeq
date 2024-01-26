
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

## Analysis of Uracils Created by human APOBEC3A and APOBEC3B*

**Step1: Download the raw sequences using the provided accession numbers.<br>**
##### Samples.xlsx :  lists of all of the samples and their accession numbers.

**Four APOBEC3A samples are already published <br>**
A3A_F6&emsp;SRR6924522&emsp;PRJNA448166<br>
A3A_H1&emsp;SRR9864913&emsp;PRJNA448166<br>
A3A_A&emsp;SRR17822878&emsp;PRJNA801888<br>
A3A_G&emsp;SRR17822877&emsp;PRJNA801888<br>
**New samples**<br>
3 A3B-CTD + 3 Empty vector control samples<br>
6 A3B-full + 6 Empty vector control samples
______
**Step2: Sequence alignment on and extracting the depth of coverage and the nucleotide counts tables.<br>**
**The reference sequences are available inside the BH214_genome_files directory.<br><br>**
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
## Hairpin Signature Analysis
**Run the Hairpin Signature Analysis in MATLAB<br>**
**The example MATLAB code is also available at Hairpin_Signature_Analysis/Hairpin_analysis.m<br>**
***To Run hairpin_signature_analysis function, the input table(struct) needs the following fields:<br>**
_ref:_ reference sequence, numeric : A=1, C=2, G=3, T=4<br>
_alt:_ alternative sequence, numeric : A=1, C=2, G=3, T=4<br>
_looplen:_ length of the potential hairpin loop<br>
_looppos:_ position of the cytosine in the potential hairpin loop<br>
_ss:_ stem strength of the potential hairpin loop<br>
_minus0:_ base immediately 5' to the cytosine, numeric : A=1, C=2, G=3, T=4<br>
_patid:_ patient/sample index, numeric (optional)<br>

***To extract the hairpin information and create a table with the above fields, you can use get_hairpin_info function.<br>**
***get_hairpin_info function requires path to the reference fasta file, and a table of mutations with the following fields:<br>**
_chr:_ chromosome<br>
_pos:_ position of the mutation<br>
_ref:_ reference sequence, numeric : A/C/G/T or 1/2/3/4<br>
_alt:_ alternative sequence, numeric : A/C/G/T or 1/2/3/4<br><br>
***samtools must be available at the system level to run get_hairpin_info function<br>** 
***Alternatively, one can extract the hairpin informations from the output of ApoHP tool<br>**
***Dependency functions are available at the Hairpin_Signature_Analysis folder<br>**

##### Try these functions on the Test_Mutations:<br>

```
## in bash:
$ samtools --version  ## make sure samtools is available
$ matlab -nodesktop

## in MATLAB:
addpath(cd /path/to/matlab_functions/Hairpin_Signature_Analysis)
% load the mutations
>> M = load_struct('./Test_Mutations.txt');
>> M.pos = str2double(M.pos); 
>> M.chr(strcmp(M.chr,'23'))={'X'}; 
>> M.chr(strcmp(M.chr,'24'))={'Y'};

% get the hairpin info for mutations
>> ref_fasta =  '/path/to/refernece/sequence/such/as/Homo_sapiens_assembly19.fasta';
>> Hairpins = get_hairpin_info(M,ref_fasta); % about 40 minutes for test mutations

% Hairpin Signature analysis
>> Hairpins.pat_id = arrayfun(@(x) str2num(x{1}), Hairpins.pat_id);
>> HSA = hairpin_signature_analysis(Hairpins,true); %% second argument for filtering only TpC sites.
>> pr(HSA.pat)
% pat_id hs1   hs2   log2R   judgement
% 1      0.046 0.037 0.33    A3A-like 
% 2      0.014 0.021 -0.56   A3B-like 
% 13     0.019 0.019 -0.0039 -        
% 14     0.020 0.021 -0.051  A3B-like 
% 50     0.061 0.043 0.51    A3A-like 

% save the output to a text file
>> save_struct(HSA.pat, 'Hairpin_Analysis.txt')

```
<br>

**Note:** <br>
For large number of mutations it might be faster to servey the entire genome for hairpins.<br>
Then map the mutations to the hairpins sites and extract looplen, looppos and ss values.<br>
See [https://github.com/alangenb/ApoHP](https://github.com/alangenb/ApoHP) for ApoHP: APOBEC hairpins analysis tool

________________
### If you have any questions regarding the article or the code behind it, please contact the corresponding author or Ramin Sakhtemani rsakhtemani (at) mgh.harvard.edu
