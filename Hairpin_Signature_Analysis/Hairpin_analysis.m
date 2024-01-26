%% bash 
%% Make sure samtools is available on the system
cd path/to/Hairpin_Signature_Analysis
module load matlab/2016b
matlab -nodesktop
%%%%%%%%%%%%%%%%%%%%%
%% matlab:
% load the mutations file
% replace the path/to/Test_Mutations.txt with the path to the mutations file
M = load_struct('path/to/Test_Mutations.txt');
M.pos = str2double(M.pos); 
M.chr(strcmp(M.chr,'23'))={'X'}; 
M.chr(strcmp(M.chr,'24'))={'Y'};

%%%%%%%%%%%%%%%%%%%%%
%% ref_fasta is the path to the reference fasta file.
%% make sure to use the correct fasta file.
%% Test_Mutations.txt is based on hg19.
ref_fasta =  '/path/to/reference.fa/Homo_sapiens_assembly19.fasta';
%%%%%%%%%%%%%%%%%%%%%
%% get the hairpin info for mutations 
Hairpins = get_hairpin_info(M,ref_fasta); % about 40 minutes
%%%%%%%%%%%%%%%%%%%%%
%% Hairpin Signature analysis
%% note: patient/sample id column must be numeric
Hairpins.pat_id = arrayfun(@(x) str2num(x{1}), Hairpins.pat_id);
HSA = hairpin_signature_analysis(Hairpins,true); %% second argument for filtering only TpC sites.
% >> pr(HSA.pat)
% pat_id hs1   hs2   log2R   judgment
% 1      0.046 0.037 0.33    A3A-like 
% 2      0.014 0.021 -0.56   A3B-like 
% 13     0.019 0.019 -0.0039 -        
% 14     0.020 0.021 -0.051  A3B-like 
% 50     0.061 0.043 0.51    A3A-like 
%% save the output to a text file
save_struct(HSA.pat, 'Hairpin_Analysis.txt')

