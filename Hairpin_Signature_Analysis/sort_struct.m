function [s2,ord]=sort_struct(s1,keyfield,order)
%
% sort_struct(struct,keyfield,order)
%
% sorts structure by specified keyfield
%    in descending order if order=-1
%
% keyfield and order can be arrays, for nested sorting
%
% Mike Lawrence 2008-05-01

if slength(s1)==0
  s2 = s1;
  ord = [];
  return
end

if length(keyfield)==0, return; end

if ~iscell(keyfield)
  keyfield = {keyfield};
end

if ~exist('order','var')
  order = repmat(1,length(keyfield),1);
end

if ischar(order) && strcmpi(order,'descend')
  order = [-1];
end

if length(order) ~= length(keyfield)
  error('order and keyfield must have same number of elements');
end

if any(order~=1 & order~=-1) error('unknown order type'); end

orig_len = slength(s1);
ord=(1:orig_len)';
fields = fieldnames(s1);
nf = length(fields);

rank = nan(orig_len,nf);

for k=1:length(keyfield)
  f = getfield(s1,keyfield{k});
  if length(f)<orig_len, error('Attempted to sort on truncated field "%s"',keyfield{k}); end
  if islogical(f), f=1*f; end
  if isnumeric(f)
    [u ui uj] = unique(f,'rows');
    [tmp ordi] = sortrows(u);   
  else
    [u ui uj] = unique(f);
    [tmp ordi] = sort(u);
  end
  if order(k)==-1, ordi=ordi(end:-1:1); end
  rank(:,k) = ordi(uj);
end

[tmp ord] = sortrows(rank);

s2 = reorder_struct(s1,ord);
