function x = get_first_lines(fname,n)

if ~exist('n','var'), n=1; end

x = cell(n,1);

try
  f = fopen(fname,'rt');
  for i=1:n
    tmp = fgetl(f);
    if tmp==-1, break; end
    x{i} = tmp;
  end
  fclose(f);
catch me
  fclose(f);
end



