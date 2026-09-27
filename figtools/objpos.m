function [xy] = objpos(obj,pp,opts)
arguments
    obj
    pp
    opts.axesPlotBox = true
end

switch(obj.Type)
    case 'text'
        pos = obj.Extent;
    case 'figure'
        pos = obj.Position;
    case 'axes'
        switch(opts.axesPlotBox)
            case true
                pos = getAxesTruePosition(obj);
            case false
                pos = obj.Position;
        end
    otherwise
        pos = [min(obj.XData(:)) min(obj.YData(:))];
        pos = [pos [max(obj.XData(:)) max(obj.YData(:))]-pos];
end

cent = @(x)x*[1 0 0.5 0]';
mid = @(x)x*[0 1 0 0.5]';

switch(pp)
    case {'bot','b','s'}
        xy = [cent(pos) pos(2)];

    case {'top','t','n'}
        xy = [cent(pos) pos*[0 1 0 1]'];

    case {'left','l','w'}
        xy = [pos(1) mid(pos)];

    case {'right','r','e'}
        xy = [pos*[1 0 1 0]' mid(pos)];

    case {'center','c'}
        xy = [cent(pos) mid(pos)];

    case {'botleft','bl','sw'}
        xy = [pos(1:2)];

    case {'topleft','tl','nw'}
        xy = [pos(1) pos*[0 1 0 1]'];

    case {'botright','br','se'}
        xy = [pos*[1 0 1 0]' pos(2)];

    case {'topright','tr','ne'}
        xy = [pos*[1 0 1 0]' pos*[0 1 0 1]'];

    otherwise
        xy = [pos*[1 0 pp(1) 0]' pos*[0 1 0 pp(end)]'];

end


end