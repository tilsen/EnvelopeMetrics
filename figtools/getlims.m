function [lims] = getlims(ax,xy)

if nargin==1
    xy = 'xy';
end

switch(xy)
    case 'y'
        lims = get(ax,'ylim'); 
        if length(ax)>1, lims = minmax(reshape(vertcat(lims{:}),1,[])); end
    case 'x'
         lims = get(ax,'xlim');
        if length(ax)>1,  lims = minmax(reshape(vertcat(lims{:}),1,[])); end
    case 'xy'
        if length(ax)>1
            lims = get(ax,'ylim');  limsy = minmax(reshape(vertcat(lims{:}),1,[]));
            lims = get(ax,'xlim');  limsx = minmax(reshape(vertcat(lims{:}),1,[]));
            lims = [limsx; limsy];
        else
            lims = [get(ax,'xlim'); get(ax,'ylim')];
        end
        
end

lims = double(lims);

end