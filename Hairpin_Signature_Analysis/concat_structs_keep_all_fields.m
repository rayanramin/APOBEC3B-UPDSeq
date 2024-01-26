function X = concat_structs_keep_all_fields(X)
% Mike Lawrence 2010

% get comprehensive list of fields
allflds = [];
type = [];  % 1 = numeric, 2 = boolean, 3 = cell
szs = {};  % dimensions (first dimension zeroed out)
for i=1:length(X)
  fn = fieldnames(X{i});
  allflds = [allflds; fn];
  for j=1:length(fn)
    g = getfield(X{i},fn{j});
    sz = size(g); sz(1)=0; sz = sprintf('%d ',sz);
    szs = [szs; sz];
    if isnumeric(g), t = 1;
    elseif islogical(g), t=2;
    elseif iscell(g), t = 3;
    else error('Unknown type: input %d field %s\n',i,fn{j});
    end
    type = [type; t];
  end
end
[flds ui uj] = unique(allflds,'first');
%count(allflds);

% make sure types are compatible
ot = type;
os = szs;
type = nan(length(flds),1);
szs = cell(length(flds),1);
for i=1:length(flds)
  t = ot(uj==i);
%  if length(unique(type(i)))>1, error('Inputs have "%s" of different types',flds{i}); end    % (this error check has never been correctly imposed)
  type(i) = t(1);
  s = os(uj==i);
  if length(unique(s))>1, error('Inputs have "%s" of different dimensions',flds{i}); end
  szs(i) = s(1);
end

% for each input, if it doesn't have the field in question, then add a blank version
for i=1:length(X)
  for j=1:length(flds)
    if ~isfield(X{i},flds{j})
      sz = as_row(sscanf(szs{j},'%d ')); sz(1)=slength(X{i});
      if type(j)==1, z = nan(sz);
      elseif type(j)==2, z = false(sz);
      elseif type(j)==3, z = repmat({''},sz);
      else error('Inconsistent behavior');
      end
      X{i} = setfield(X{i},flds{j},z);
    end
  end
end

X = concat_structs(X);
