function truePos = getAxesTruePosition(ax)
% getAxesTruePosition  Returns the true rendered position of an axes,
%   accounting for PlotBoxAspectRatio, DataAspectRatio (axis equal /
%   axis image), or neither. Correctly handles non-square figures.
%
% Usage:
%   truePos = getAxesTruePosition(ax)
%
% Input:
%   ax      - Handle to a valid axes object
%
% Output:
%   truePos - 1x4 vector [x, y, width, height] in normalised figure units

    if nargin < 1 || ~isgraphics(ax, 'axes')
        error('Input must be a valid axes handle.');
    end

    darConstrained  = strcmp(ax.DataAspectRatioMode,    'manual');
    pbarConstrained = strcmp(ax.PlotBoxAspectRatioMode, 'manual');

    if ~darConstrained && ~pbarConstrained
        truePos = ax.Position;
        return;
    end

    % Determine the container rectangle within which MATLAB centres the
    % constrained plot box.
    if isprop(ax, 'PositionConstraint')
        posConstraint = ax.PositionConstraint;     % R2020a+
    else
        posConstraint = ax.ActivePositionProperty; % pre-R2020a
    end

    if strcmp(posConstraint, 'outerposition')
        % Outer rectangle is fixed; subtract label/tick space to get the
        % inner container.
        op  = ax.OuterPosition;
        ti  = ax.TightInset;       % [left bottom right top]
        pos = [op(1) + ti(1), ...
               op(2) + ti(2), ...
               op(3) - ti(1) - ti(3), ...
               op(4) - ti(2) - ti(4)];
    else
        % 'innerposition' / 'position': user-set, use directly.
        pos = ax.Position;
    end

    % Figure pixel dimensions are required to convert between normalised
    % units and the display aspect ratio used by MATLAB's renderer.
    % A normalised width of 1 spans figSz(1) pixels; height spans figSz(2).
    fig      = ancestor(ax, 'figure');
    prevU    = fig.Units;
    fig.Units = 'pixels';
    figSz    = fig.Position(3:4);  % [width height] in pixels
    fig.Units = prevU;
    figAR    = figSz(1) / figSz(2);

    % Displayed (pixel) aspect ratio of the data / plot box.
    if pbarConstrained
        % PlotBoxAspectRatio takes precedence when both modes are manual.
        pbr         = ax.PlotBoxAspectRatio;
        aspectRatio = pbr(1) / pbr(2);
    else
        % DataAspectRatio (axis equal, axis image, …).
        % abs() guards against reversed axes (YDir='reverse' from imshow).
        dar         = ax.DataAspectRatio;
        xSpan       = abs(diff(ax.XLim)) / dar(1);
        ySpan       = abs(diff(ax.YLim)) / dar(2);
        aspectRatio = xSpan / ySpan;
    end

    % Fit the constrained plot box inside `pos` in pixel space so that
    % non-square figures are handled correctly.
    nominalAR = (pos(3) / pos(4)) * figAR;  % displayed AR of the container

    if nominalAR > aspectRatio
        % Container is wider than needed — height is the binding limit.
        height = pos(4);
        width  = pos(4) * aspectRatio / figAR;
    else
        % Container is taller than needed — width is the binding limit.
        width  = pos(3);
        height = pos(3) * figAR / aspectRatio;
    end

    xCenter = pos(1) + pos(3) / 2;
    yCenter = pos(2) + pos(4) / 2;
    truePos = [xCenter - width/2, yCenter - height/2, width, height];
end