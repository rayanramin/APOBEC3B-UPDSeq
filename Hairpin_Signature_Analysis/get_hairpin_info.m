function OUT = get_hairpin_info(M , ref_fasta)
	%%% check if M is a struct 
	if ~isstruct(M)
		error('M must be a struct')
	end
	%%% check if M has the right fields (chr, pos, ref, alt)
	if ~all(isfield(M,{'chr','pos','ref','alt'}))
		error('M must have the fields chr, pos, ref, alt' )
	end
	%%% check if M.chr is a cell array of strings
	if ~iscellstr(M.chr)
		error('M.chr must be a cell array of strings')
	end
	%%% check if M.pos is a vector of integers
	if ~isnumeric(M.pos) || ~all(M.pos==floor(M.pos))
		error('M.pos must be a vector of integers')
	end
	%%% check if ref_fasta is a string ending in .fa or .fasta
	if ~ischar(ref_fasta) || ~ grepm('\.fasta$|\.fa$|\.fna$',{ref_fasta})
		error('ref_fasta must be a string ending in .fa or .fasta')
	end
	%%% check if ref_fasta exists
	if ~exist(ref_fasta,'file') 
		error('ref_fasta does not exist')
	end
	%%% check if samtools is available on the system
	[status, ~] = system('samtools --version');
	if status ~= 0; error('samtools not available on the system');	end
	%%% check if fasta.fai exists and if not create it
	if ~exist([ref_fasta,'.fai'],'file')
		[status, ~] = system(['samtools faidx ',ref_fasta]);
		if status ==0 & exist([ref_fasta,'.fai'],'file')
			fprintf('created %s.fai\n',ref_fasta)
		else
			error('samtools faidx failed')
		end
	end
	%%% read the first column of the fai file and make sure all M.chr are in it
	[~, chr_names] = system(['cut -f1 ',ref_fasta,'.fai']);
	chr_names = strsplit(chr_names,'\n');
	m_chr = unique(M.chr);
	if ~all(ismember(m_chr, chr_names));
		error(['M.chr contains chromosomes not in the reference fasta (',chr_names{1},', ',chr_names{2},', ',chr_names{3} ,'...)'])
	end
	%%% check if the first mutation can be loaded
	st_pos = M.pos(1)-200; end_pos = M.pos(1)+200;
	[status, seq] = system(['samtools faidx ',ref_fasta,' ',M.chr{1},':',num2str(st_pos),'-',num2str(end_pos)]); 
	if status ~= 0; error('Failed check: (fetching first mutation context)');	end
	%%% check if M.ref is a vector of integers between 1 and 4 or a vector of strings A,C,G,T
	if ~isnumeric(M.ref); if all(ismember(M.ref,{'A','C','G','T'})); 
		M.ref = uint8(mapacross(upper(M.ref),{'A','C','G','T'},1:4));
	end;end
	if isnumeric(M.ref) & ~all(M.ref>=1 & M.ref<=4) 
		error('M.ref must be a vector of integers between 1 and 4 or a vector of strings A,C,G,T')
	end
	if ~isnumeric(M.alt); if all(ismember(M.alt,{'A','C','G','T'})); 
		M.alt = uint8(mapacross(upper(M.alt),{'A','C','G','T'},1:4));
	end;end
	if isnumeric(M.alt) & ~all(M.alt>=1 & M.alt<=4) 
		error('M.alt must be a vector of integers between 1 and 4 or a vector of strings A,C,G,T')
	end
	%%%%%%%%%%%%%%%%%%%%%
	%%%% Survey the mutations for hairpin information
	fprintf('Analysing genomic positions...\n');
	minlooplen=3;      % consider loops of at least 3nt
	maxlooplen=11;     % consider loops of at most 11nt
	maxstem=100;       % consider stems of up to 100bp
	minbulgepos=2;     % consider hairpins with bulges at least 2bp away from loop
	maxbulgepos=8;     % consider hairpins with bulges up to 8bp away from loop
	minmismatchpos=2;  % consider hairpins with mismatches in the stem at least 2 bp away from loop
	maxintrastem=1;    % consider positions to be intrastem if they are within 1 bp of loop
	flank = 200;
	ns = flank*2 +1;
	mask = 30 ; % to avoide changing the code
	fs = {'looplen','looppos','bulgepos','nbp','ngc','mmp','ss'};

	OUT=struct();
	revstr='';
	for mut=1:slength(M)
		if M.ref(mut) == 1| M.ref(mut)==4; continue; end % skip A/T's
		chr = M.chr{mut}; pos = M.pos(mut);
		st_pos = pos-flank; end_pos = pos+flank;
		[status, seq] = system(['samtools faidx ',ref_fasta,' ',chr,':',num2str(st_pos),'-',num2str(end_pos)]); 
		if status ~= 0; error(['Failed to get sequence for mutation #',num2str(mut)]);end
		seq = strsplit(seq,'\n');
		ref = uint8(listmap( regexprep(upper(strjoin(seq(find(grepm('^>.*',seq))+1:end),'')),'[^ACGTN]','' )','ACGT'));
		if length(ref) ~= ns; error('Something went wrong: ref length'); end
		if ref(flank+1) ~= M.ref(mut); error('Something went wrong: ref identity'); end
		XSite= struct();
		XSite.ref = ref ;
		XSite.minus2  = [0;0;0;XSite.ref(1:end-3)];
		XSite.minus1  = [0;0;XSite.ref(1:end-2)];
		XSite.minus0  = [0;XSite.ref(1:end-1)];
		XSite.plus1   = [XSite.ref(2:end);0];
		XSite.plus2   = [XSite.ref(3:end);0;0];
		XSite.plus3   = [XSite.ref(4:end);0;0;0];
		for i=1:length(fs), XSite.(fs{i}) = -ones(ns,1,'int8'); end
		gc = int8(XSite.ref==3 | XSite.ref==2); % (makes things quicker in loop)
		
		% fprintf('Loop length: ');
		for looplen=minlooplen:maxlooplen, % fprintf(' %d/%d',looplen, maxlooplen);
		for looppos=1-maxintrastem:looplen+maxintrastem
			for bulgepos=[0 (-minbulgepos):-1:-(maxbulgepos) minbulgepos:maxbulgepos]  % 0 = no bulge
			stem_open = false(ns,1); stem_open((1+mask):(ns-mask))=true; can_pair = false(ns,1);
			nbp = zeros(ns,1,'int8'); ngc = zeros(ns,1,'int8'); mmp = zeros(ns,1,'int8');
			lft=looppos; rgt=1+looplen-looppos; done=false;
			k_switch_to_list_mode = 3;  % after first two basepairs, switch from global mode to list mode (more efficient)
			for k=1:maxstem
				if k<k_switch_to_list_mode    % optimized for early steps (when all/most of genome is under consideration)
					can_pair((1+mask):(ns-mask)) = XSite.ref((1+mask-lft):(ns-mask-lft)) == 5-XSite.ref((1+mask+rgt):(ns-mask+rgt));
					extend = (stem_open & can_pair);
					nbp(extend)=nbp(extend)+1;
					tmp=[zeros(lft,1,'int8');gc(1:end-lft)]; ngc(extend) = ngc(extend) + tmp(extend);
					if k>minmismatchpos
						stem_open(~can_pair & mmp>0)=0;
						extend_over_mismatch = (stem_open & ~can_pair & mmp==0);
						mmp(extend_over_mismatch)=k;
					else
						stem_open(~can_pair)=0;
					end
				done = ~any(stem_open);
				else      % same procedure optimized for later steps (when smaller lists of positions remain)
					ii = find(stem_open);
					can_pair = XSite.ref(ii-lft) == 5-XSite.ref(ii+rgt);
					fe = ii(can_pair);
					if ~isempty(fe)
						nbp(fe)=nbp(fe)+1;
						ngc(fe)=ngc(fe)+gc(fe-lft);
					end
					if k>minmismatchpos
						stem_open(ii(~can_pair & mmp(ii)>0))=0;
						feomm = ii(~can_pair & mmp(ii)==0);
						mmp(feomm) = k;
					else
						stem_open(ii(~can_pair))=0;
					end
				done = ~any(stem_open(ii));
				end
				lft=lft+1; if k==-bulgepos, lft=lft+1; end            % bypass left bulge
				rgt=rgt+1; if k==+bulgepos, rgt=rgt+1; end            % bypass right bulge
				if done,break;end
			end
			
			% remove terminal mismatches
			mmp(mmp==nbp+1) = 0;
			% score hairpin strength
			ss = nbp + 2*ngc;                                       % base strength
			ss = ss + max(-100*mmp, -max(1,13-2*mmp));              % mismatch penalty
			if bulgepos~=0, ss = ss - 6; end                        % bulge penalty
			% assign hairpins, overwriting weaker with stronger
			idx = find(ss>XSite.ss);
			for i=1:length(fs), f=fs{i}; tmp=eval(f); if length(tmp)==1, XSite.(f)(idx) = tmp; else XSite.(f)(idx) = tmp(idx); end; end
			end
		end
		end
		% for G/A positions on reference strand, flip annotations to opposite strand so that we're always annotating C/T's
		ga = (XSite.ref==3 | XSite.ref==1);
		tmp=XSite.plus1(ga); XSite.plus1(ga)=5-XSite.minus0(ga); XSite.minus0(ga)=5-tmp;
		tmp=XSite.plus2(ga); XSite.plus2(ga)=5-XSite.minus1(ga); XSite.minus1(ga)=5-tmp;
		tmp=XSite.plus3(ga); XSite.plus3(ga)=5-XSite.minus2(ga); XSite.minus2(ga)=5-tmp;
		XSite.looppos(ga) = XSite.looplen(ga)-XSite.looppos(ga)+1;
		XSite.bulgepos(ga) = -XSite.bulgepos(ga);
		
		tmp = struct();
		tmp.chr = {chr}; tmp.pos = pos; 
		tmp.alt = M.alt(mut);
		%%% if mut_id and pat_id or samp_id are present in M, add them to tmp
		fn=fieldnames(M);
		if any(ismember(fn,{'pat_idx','patid','patidx','patient_id','sample_id','samp_id','samp_idx','IDX'}));
			tmp.pat_id = M.(fn{ismember(fn,{'pat_idx','patid','patidx','patient_id','sample_id','samp_id','samp_idx','IDX'})})(mut);
		end
		if any(ismember(fn,{'mut_idx','mutid','mutidx','mutation_id','mutation_idx','mutIDX'}));
			tmp.mut_id = M.(fn{ismember(fn,{'mut_idx','mutid','mutidx','mutation_id','mutation_idx','mutIDX'})})(mut);
		end
		if any(ismember(fn,{'samp_idx','sampid','sampidx','sample_id','sample_idx','sampleIDX'}));
			tmp.samp_id = M.(fn{ismember(fn,{'samp_idx','sampid','sampidx','sample_id','sample_idx','sampleIDX'})})(mut);
		end
		%%% add the Hairpin info to tmp
		tmp=merge_structs({tmp,reorder_struct(XSite,flank+1)});

		if mut==1
			OUT=tmp;
		else 
			OUT=concat_structs_keep_all_fields({OUT, tmp}); 
		end

		pct = 100*mut/slength(M);
		msg = sprintf('Analyzing Hairpins: %3.1f',pct);
		fprintf([revstr, msg]);
		revstr = repmat(sprintf('\b'),1,length(msg));
	end

	OUT = move_field_to_before(OUT,'ref','pat_id');

	% %%%% check to see if the length of OUT is the same as the length of M
	if slength(OUT) == sum(M.ref==3 | M.ref==2)
		fprintf('\nFinal check : OK\n')
	else
		warning('\nSomething went wrong: OUT length\n')
	end

