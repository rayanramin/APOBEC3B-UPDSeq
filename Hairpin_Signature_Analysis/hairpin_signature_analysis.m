function OUT =  hairpin_signature_analysis(X, only_TpC)
	%%% if only_TpC is not provided, default value for only_TpC is false
	if ~exist('only_TpC', 'var')
		only_TpC = true
	end

	%%% make sure input X is provided
	if nargin < 1
		error('Not enough input arguments')
	end
	if ~isstruct(X)
		error('input must be a struct')
	end
	%%% check if M has the right fields  
	if ~all(isfield(X,{'looplen','looppos','ss','minus0'}))
		error('input struct must have the fields looplen, looppos, ss,minus0' )
	end
	%%% 
	if only_TpC & ~all(isfield(X,{'minus0'}))
		error('input struct must have the field minus0' )
	end
	%%%%% 
	%%% bin the ss values {'0-4','5-7','8-11','12-15','16+'}
	fprintf('binning ss values')
	[~, X.bin]=histc(X.ss ,[0;5;8;12;16;Inf]) ;
	fprintf([repmat(sprintf('\b'),1, 17),'binning ss values: done\n' ])
	%%%%%%% check if X has 'pat_id' or 'patid' or 'patient_id' fields
	if sum(isfield(X,{'pat_id','pat_idx','patid','patidx','patient_id','sample_id','samp_id','samp_idx','IDX'}))>1
		error('X must have only one of the fields pat_id, pat_idx, patid, patient_id, IDX')
	else 
		if sum(isfield(X,{'pat_id','pat_idx','patid','patidx','patient_id','sample_id','samp_id','samp_idx','IDX'}))==0
			X.pat_id=repmat(1,slength(X),1);
		end
	end
	%%% get the list of pat_idx
	%%% if pat_idx is not provided then rename the other fields to pat_idx
	if ~isfield(X,'pat_id')
		fn=fieldnames(X);
		X.pat_id = X.(fn{ismember(fn,{'pat_idx','patid','patidx','patient_id','sample_id','samp_id','samp_idx','IDX'})});
	end
	if ~isnumeric(X.pat_id)
		error('pat_id must be numeric')
	end
	
	%%%% %%%% %%%%
	OUT = struct();
	OUT.pat = struct();
	[OUT.pat.pat_id, ~, ~] = unique(X.pat_id);
	%%%
	
	c2gt = (X.ref==2 & X.alt~=1)|(X.ref==3 & X.alt ~=4); % C>T or G>A
	nC2GT = histc(X.pat_id(c2gt),OUT.pat.pat_id);
	tc2gt = (c2gt & X.minus0==4 ); % TpC
	nTC2GT = histc(X.pat_id(tc2gt),OUT.pat.pat_id);
	
	if only_TpC & any(isfield(X,{'minus0'}))
		X = reorder_struct(X,X.minus0==4);
		fprintf('filtered for TpC positions only\n')
	end
	
	
	%%% let's count the mutations
	fprintf('counting hairpins')
	C = zeros(length(OUT.pat.pat_id), 90 );
	looplens = [3 4 5 6];
 	X.HPGroup = zeros(slength(X),1);

	HP = struct();
	HP.hairpin_group = [1:90]';
	HP.looplen = [repmat(3,3*5,1);repmat(4,4*5,1);repmat(5,5*5,1);repmat(6,6*5,1)];
	HP.looppos = [repelem([1:3]',5);repelem([1:4]',5);repelem([1:5]',5);repelem([1:6]',5)];
	HP.ssbin = repmat([1:5]',18,1);
	HP.ss_min = mapacross(HP.ssbin,[1:5]',[0;5;8;12;16]);
	HP.ss_max = mapacross(HP.ssbin,[1:5]',[4;7;11;15;Inf]);
	
	HP.HS1=[13,15,5,0,0,37,37,18,9,2,64,131,175,246,111,...
	38,23,13,3,1,35,33,13,6,1,37,60,28,20,32,...
	30,81,69,79,119,14,16,8,3,0,29,36,21,1,1,...
	25,37,19,8,0,21,50,48,49,72,26,41,24,3,12,...
	21,12,9,0,1,11,13,9,3,3,22,47,12,7,1,...
	19,43,16,7,12,21,49,33,17,39,31,36,24,3,2]' .* 0.8; % 0.8 scale factor to make nmf work better
	HP.HS2=[59,66,24,9,1,82,89,42,28,2,194,161,121,128,74,...
	66,85,41,3,1,101,96,41,8,3,114,126,43,20,8,...
	105,222,524,549,344,34,62,28,7,3,69,89,40,9,8,...
	81,83,57,21,4,62,116,73,49,37,95,162,159,176,119,...
	44,55,23,2,1,62,73,41,14,2,40,108,62,14,15,...
	34,98,55,26,8,46,119,91,47,44,71,122,72,30,36]';
	
	X.HPGroup = zeros(slength(X),1);
	for i=1:90
		hp_use = X.looplen == HP.looplen(i) & X.looppos == HP.looppos(i) & X.bin == HP.ssbin(i);
		if sum(hp_use) ~= 0;
			X.HPGroup(hp_use) = HP.hairpin_group(i);
		end
		if i==90; fprintf([repmat(sprintf('\b'),1, 17),'counting hairpins: done\n' ]); end
	end
	

	%%% count per patients
	revstr = '';
	for p = 1:slength(OUT.pat);
		C(p,:) = histc(X.HPGroup(X.pat_id == OUT.pat.pat_id(p)),[1:90]);
		pct = 100*p/slength(OUT.pat);
		msg = sprintf('counting per patient: %3.1f',pct);
		fprintf([revstr, msg]);
		revstr = repmat(sprintf('\b'),1,length(msg));
	end
	fprintf([repmat(sprintf('\b'),1, 5),' done\n' ]);
	
	fprintf('running NMF')
	damp=0.03;
	hs = [HP.HS1' ; HP.HS2'];
	N = nmf3(C+1,hs);
	r = log2((N(:,1)+damp)./(N(:,2)+damp));
	fprintf([repmat(sprintf('\b'),1, 11),'running NMF: done\n' ]);

	fprintf('assigning signatures')
	threshold = 0.02; min_nC =200; min_tc_frac = 0.20259*1.75; % ~50% more TpC than random;


	tc_frac = nTC2GT./nC2GT ;
	OUT.pat.nC2GT = nC2GT; OUT.pat.nTC2GT = nTC2GT; OUT.pat.tc_frac = tc_frac;
	
	OUT.pat.hs1 = N(:,1);
	OUT.pat.hs2 = N(:,2);
	OUT.pat.log2R = r;
	OUT.pat.judgment=repmat({'-'} ,length(OUT.pat.pat_id),1);
	OUT.pat.judgment(r>threshold & nC2GT >= min_nC & tc_frac >= min_tc_frac) = {'A3A-like'};
	OUT.pat.judgment(r< -threshold & nC2GT >= min_nC & tc_frac >= min_tc_frac) = {'A3B-like'};

	OUT.counts = C;
	OUT.HP = HP;
	if only_TpC
		OUT.only_TpC = true;
	else
		OUT.only_TpC = false;
	end

	fprintf(': done!\n')
