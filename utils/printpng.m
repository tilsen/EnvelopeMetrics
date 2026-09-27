function [] = printpng(str,varargin)

p = inputParser;

addRequired(p,'str',@(x)ischar(x) | (isstring(x) & isscalar(x)));
addParameter(p,'resolution','-r400');

str = char(str);
parse(p,str,varargin{:});

res =p.Results;

print(gcf,'-dpng',res.resolution,[str '.png']);

end