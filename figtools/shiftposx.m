function [] = shiftposx(obj,shift)

if isscalar(shift)
    shift = repmat(shift,1,numel(obj));
end

for i=1:numel(obj)
    if ~ishandle(obj(i)),continue; end
    switch(obj(i).Type)
        case {'patch' 'line'}
            obj(i).XData = obj(i).XData + shift(i);
        otherwise
            obj(i).Position(1) =  obj(i).Position(1) + shift(i);
    end
end
end