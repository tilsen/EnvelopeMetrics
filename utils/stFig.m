classdef stFig < handle

    properties
        Handle
        Parent = [] % explicit attachment; empty retains standalone construction
        Axes
        PanelSpec = [1 1]
        TensorSpec
        PanelLayout
        Panels
        ExteriorMargins = [0.05 0.05 0.01 0.01] %left, bottom, right, top
        InteriorMargins = [0.05 0.05]
        AxisPositions
        AxesLabels = gobjects(0);
        Toolbar = 'none'
        Menubar = 'none'
        Name = ""
        NumberTitle = "off"
        Scrollable = "off"
        Units = "Normalized"
        Position = [] % [] resolved at figure creation: primary-monitor-size OuterPosition if running in a Live Script, MATLAB's own default placement otherwise (see createFigure)
        Aspect = "Maximized"
        AxesHandleForm = 'vector'        
        Tags = []        
        Theme = "light"
        MakeBackgroundAxes = false;
        BackgroundAxes
        OverlayAxes = gobjects(0) % creation order; newest live overlay is the default
        AxesTickProps = {'TickDir','out','TickLen',0.003*[1 1]};
        UserData = {}
        CurrentAxesIx = []
        currentAxes = [];
    end

    methods
        function obj = stFig(panelSpec,ExteriorMargins,InteriorMargins,opts)
            arguments
                panelSpec = [1 1]
                ExteriorMargins = []
                InteriorMargins = []
                opts.Parent = []
                opts.Theme = []
                opts.Position = []
                opts.Aspect = []
                opts.AxesHandleForm = []
                opts.NumberTitle = [];
                opts.Scrollable = [];
                opts.Name = [];
                opts.MakeBackgroundAxes = []
                opts.TensorSpec = []
                opts.Tags = []
            end

            obj.PanelSpec = panelSpec;
            if ~isempty(ExteriorMargins), obj.ExteriorMargins = ExteriorMargins; end
            if ~isempty(InteriorMargins), obj.InteriorMargins = InteriorMargins; end

            %handle simplified margin inputs
            if numel(obj.ExteriorMargins)==1 %#ok<ISCL>
                obj.ExteriorMargins = obj.ExteriorMargins*[1 1 1 1];
            elseif numel(obj.ExteriorMargins)==2
                obj.ExteriorMargins = repmat(obj.ExteriorMargins(:)',1,2);
            end            
            if numel(obj.InteriorMargins)==1 %#ok<ISCL>
                obj.InteriorMargins = obj.InteriorMargins*[1 1];
            end

            ff = fieldnames(opts);
            for i=1:length(ff)
                if ~isempty(opts.(ff{i}))
                    obj.(ff{i}) = opts.(ff{i});
                end
            end

            switch(obj.AxesHandleForm)
                case 'tensor'
                    [obj.PanelLayout,obj.Panels] = obj.processPanelSpec(obj.PanelSpec(1:2));
                case 'transposed'
                    [obj.PanelLayout,obj.Panels] = obj.processPanelSpec(obj.PanelSpec);
                    obj.PanelLayout = obj.PanelLayout';
                    obj.Panels = obj.Panels';
                otherwise
                    [obj.PanelLayout,obj.Panels] = obj.processPanelSpec(obj.PanelSpec);

            end

            if isempty(obj.Parent)
                obj.createFigure();
            else
                if isa(obj.Parent,'stFig'), obj.Parent = obj.Parent.Handle; end
                if all(isgraphics(obj.Parent,'axes'))
                    obj.Handle = ancestor(obj.Parent(1),'figure');
                    if any(arrayfun(@(a)ancestor(a,'figure')~=obj.Handle,obj.Parent))
                        error('stFig:Parent','Parent axes must belong to one figure.');
                    end
                elseif isscalar(obj.Parent) && isgraphics(obj.Parent,'figure')
                    obj.Handle = obj.Parent;
                    if ~isempty(obj.Handle.CurrentAxes)
                        obj.Parent = obj.Handle.CurrentAxes;
                    end
                else
                    error('stFig:Parent','Parent must be live axes or a figure.');
                end
            end
            obj.createAxes();

            switch(obj.AxesHandleForm)
                case {'matrix','array','grid','tensor'}
                    form = unique(obj.PanelLayout,"rows");
                    form = unique(form',"rows")';
                    ax = gobjects(size(form));

                    n=1;
                    for r=1:size(ax,1)
                        for c=1:size(ax,2)
                            ax(r,c) = obj.Axes(n);
                            n=n+1;
                            if n>numel(obj.Axes), break; end
                        end
                        if n>numel(obj.Axes),break; end
                    end
                    obj.Axes = ax;
            end

            switch(obj.AxesHandleForm)
                case 'tensor'
                    for i=1:numel(obj.Axes)
                        obj.subAxes(i,obj.TensorSpec{:});
                    end
            end

            if isempty(opts.Tags)
                obj.Tags = repmat("",size(obj.Axes));
            elseif isscalar(opts.Tags)
                obj.Tags = repmat(opts.Tags,size(obj.Axes));
            else
                obj.Tags = reshape(opts.Tags(1:numel(obj.Axes)),size(obj.Axes));
            end

            axes(obj.Axes(1));
            obj.CurrentAxesIx = 1;
            obj.currentAxes = obj.Axes(1);

        end

        function [] = createAxes(obj)
            layoutProps = obj.layoutProperties(obj.PanelLayout);

            obj.AxisPositions = obj.panelLayoutToCoords(obj.PanelLayout,obj.ExteriorMargins,obj.InteriorMargins);
            existing = gobjects(0);
            container = obj.Handle;
            if ~isempty(obj.Parent) && all(isgraphics(obj.Parent,'axes'))
                existing = obj.Parent(:);
                if numel(existing)>1 && numel(existing)~=size(obj.AxisPositions,1)
                    error('stFig:AxesCount','Supply one axes or one axes per panel.');
                end
                if isscalar(existing)
                    container = existing.Parent;
                    oldUnits = existing.Units;
                    existing.Units = 'normalized';
                    region = existing.Position;
                    existing.Units = oldUnits;
                    % The supplied panel already includes its exterior margins.
                    obj.AxisPositions = obj.panelLayoutToCoords(obj.PanelLayout,[0 0 0 0],obj.InteriorMargins);
                    obj.AxisPositions(:,1:2) = region(1:2)+obj.AxisPositions(:,1:2).*region(3:4);
                    obj.AxisPositions(:,3:4) = obj.AxisPositions(:,3:4).*region(3:4);
                end
            end
            obj.Axes = gobjects(size(obj.AxisPositions,1),1);

            for i=1:numel(obj.Axes)
                if i<=numel(existing)
                    obj.Axes(i) = existing(i);
                    if isscalar(existing)
                        set(obj.Axes(i),'Units','normalized','Position',obj.AxisPositions(i,:));
                    end
                else
                    obj.Axes(i) = axes(container,'Position',obj.AxisPositions(i,:), ...
                        'PositionConstraint','innerposition',obj.AxesTickProps{:});
                end
                if ~isempty(obj.Axes(i).Toolbar), obj.Axes(i).Toolbar.Visible = 'off'; end
                %obj.Axes(i) = axes('Position',obj.AxisPositions(i,:), obj.AxesTickProps{:});
                hold(obj.Axes(i),'on');
                %axis(obj.Axes(i),'manual');
                obj.Axes(i).UserData = {layoutProps(i),obj};
            end
            if obj.MakeBackgroundAxes
                obj.createBackgroundAxes;                
            end
        end

        function [] = createBackgroundAxes(obj)
            obj.BackgroundAxes = axes(obj.Handle,'Position',[0 0 1 1],'XColor','none','YColor','none','Clipping','off');
            hold(obj.BackgroundAxes,'on');
            plot(obj.BackgroundAxes,[0 1],[0 1],'.','color',0.999*[1 1 1],'Clipping','off');
            set(obj.BackgroundAxes,'XLim',[0 1],'YLim',[0 1]);
            uistack(obj.BackgroundAxes,'bottom');
        end

        function h = annotate(obj,endpoint,label,alignment,varargin)
            % Instance API: fig.annotate(endpoint,label,alignment,...).
            % Uses the newest overlay unless Overlay is explicitly supplied.
            h = stFig.DrawAnnotation(obj,endpoint,label,alignment,varargin{:});
        end

        function ax = createOverlay(obj)
            % Transparent figure-coordinate layer; does not select a plot axes.
            previousAxes = obj.Handle.CurrentAxes;
            restore = onCleanup(@()set(obj.Handle,'CurrentAxes',previousAxes)); %#ok<NASGU>
            ax = axes('Parent',obj.Handle,'Units','normalized', ...
                'Position',[0 0 1 1],'PositionConstraint','innerposition', ...
                'Color','none','Visible','off','Clipping','off', ...
                'XLim',[0 1],'YLim',[0 1],'XLimMode','manual','YLimMode','manual', ...
                'NextPlot','add','HitTest','off','PickableParts','none', ...
                'Tag','stFigOverlay');
            obj.OverlayAxes = obj.OverlayAxes(isgraphics(obj.OverlayAxes,'axes'));
            obj.OverlayAxes(end+1) = ax;
            uistack(ax,'top');
        end

        function ax = getOverlay(obj,ax)
            % Resolve the newest overlay, creating one on first use.
            if nargin<2 || isempty(ax)
                obj.OverlayAxes = obj.OverlayAxes(isgraphics(obj.OverlayAxes,'axes'));
                if isempty(obj.OverlayAxes)
                    ax = obj.createOverlay();
                else
                    ax = obj.OverlayAxes(end);
                end
            elseif ~isscalar(ax) || ~isgraphics(ax,'axes') || ...
                    ~any(obj.OverlayAxes==ax)
                error('stFig:InvalidOverlay','Overlay must be a live axes created by this stFig.createOverlay().');
            end
        end

        function [str] = handlePrintInput(~,str)
            if isempty(str)
                caller = dbstack;
                str = caller(3).name;
            else
                [p,f] = fileparts(str);
                if f==""
                    caller = dbstack;
                    f = caller(3).name;
                    str = p+filesep+f;
                end
            end
        end

        function [] = copyPrint(obj,str,opts)
            arguments
                obj
                str = []
                opts.padding = 'figure'
                opts.resolution = 350
            end
            str = obj.handlePrintInput(str);
            printpng(str,'resolution',"-r" + opts.resolution);
            copygraphics(obj.Handle,'Resolution',opts.resolution, ...
                'ContentType','image','padding',opts.padding);
        end
        function [] = copy(obj,opts)
            arguments
                obj
                opts.padding = 'figure'
                opts.resolution = 350
                opts.pause = false
            end
            copygraphics(obj.Handle,'Resolution',opts.resolution, ...
                'ContentType','image','padding',opts.padding);
            if opts.pause
                waitforbuttonpress;
            end
        end
        function [] = print(obj,str,opts)
            arguments
                obj
                str = []
                opts.resolution = '-r350';
                opts.format = 'png';
            end
            str = obj.handlePrintInput(str);
            printpng(str,'resolution',opts.resolution);
        end

        function [] = createFigure(obj)

            warning('off','MATLAB:uicontainer:ScrollableOnWithTextScaling');

            % Applies to axes created from here on, including ones MATLAB
            % auto-regenerates internally when new content is plotted into
            % an axes (which silently undoes a one-time `ax.Toolbar = []`,
            % so that approach can't be used per-axes at creation time).
            set(groot,'defaultAxesToolbarVisible','off');

            positionGiven = ~isempty(obj.Position);
            liveScriptSize = [];
            if ~positionGiven && obj.isLiveScript()
                % resolved here (not as a static property default) because
                % detection needs to run fresh at each figure's creation.
                % Live Script inline figure output is effectively captured
                % at primary-monitor size, so give it that size for
                % correct font/aspect fitting, matching what a real Live
                % Editor session would show. (Popups during automated
                % export()/testing runs are prevented at the MATLAB
                % process level -- launch with -noFigureWindows -- not by
                % hiding or repositioning figures here; see
                % isLiveScript()'s comments and exportReadme.m.)
                ss = get(0,'ScreenSize');
                liveScriptSize = ss(3:4);
            end

            commonArgs = {'ToolBar',obj.Toolbar,'MenuBar',obj.Menubar, ...
                'Scrollable',obj.Scrollable,'Name',obj.Name,'Theme',obj.Theme, ...
                'NumberTitle',obj.NumberTitle,'PaperPositionMode','auto','Clipping','on'};

            if positionGiven
                obj.Handle = figure('Units',obj.Units, ...
                    'OuterPosition',obj.Position, commonArgs{:});
            elseif ~isempty(liveScriptSize)
                obj.Handle = figure('Units','pixels', ...
                    'OuterPosition',[1 1 liveScriptSize], commonArgs{:});
                set(obj.Handle,'Units',obj.Units);
            else
                % Normal desktop/interactive use: no explicit position,
                % so let MATLAB place it however it normally would (its
                % own cascading default), rather than forcing a position.
                obj.Handle = figure('Units',obj.Units, commonArgs{:});
            end

            switch(isnumeric(obj.Aspect))
                case true
                    obj.setAspect();
                otherwise
                    if ~positionGiven && isempty(liveScriptSize)
                        % Only maximize when the figure wasn't already
                        % explicitly sized above. WindowState='maximized'
                        % asks the window manager for real on-screen
                        % display geometry, which is invalid when running
                        % headless under -noFigureWindows -- it produces
                        % an "Error: The figure size is not valid"
                        % warning during export()'s scene capture for a
                        % plain figure -- and is redundant anyway once
                        % OuterPosition already fixed the size
                        % (positionGiven / Live Script case above).
                        set(obj.Handle,'WindowState',obj.Aspect); drawnow;
                    end
                    % Maximizing can hand focus back to whatever figure
                    % was active before this one was created, which
                    % silently flips gcf/CurrentFigure to it -- so any
                    % later bare axes(...)/plot(...) call with no
                    % explicit Parent lands in the wrong figure. Force
                    % this figure back to current regardless.
                    set(0,'CurrentFigure',obj.Handle);
            end
        end

        function [] = setAspect(obj)
            set(obj.Handle,'units','inches');
            newpos = obj.Handle.OuterPosition;
            curr_aspect = newpos(3)/newpos(4);

            adjfac = obj.Aspect/curr_aspect;

            if curr_aspect > obj.Aspect %too wide, reduce width
                newpos(3) = newpos(3)*adjfac;
            elseif curr_aspect < obj.Aspect %reduce height
                newpos(4) = newpos(4)/adjfac;
            end

            set(obj.Handle,'outerposition',newpos); drawnow;
            set(obj.Handle,'units','normalized');

            if obj.Handle.OuterPosition(3) < 0.99
                obj.Handle.OuterPosition(1) = obj.Handle.OuterPosition(1) + (1-obj.Handle.OuterPosition(3))/2;
            end
        end

        function [] = matchAxisLims(obj,ax,dim)
            arguments
                obj
                ax = []
                dim = ''
            end

            if isempty(ax), ax = obj.Axes; end
            switch(dim)
                case 'x'
                    set(ax,'xlim',getlims(ax,'x'));
                case 'y'
                    set(ax,'ylim',getlims(ax,'y'));
                otherwise
                    set(ax,'xlim',getlims(ax,'x'));
                    set(ax,'ylim',getlims(ax,'y'));
            end

        end

        function [] = formatAxesArray(obj,opts)
            arguments
                obj
                opts.YTickLabels = 'left'
                opts.XTickLabels = 'bottom'
                opts.diagonal = ''
            end

            validAxes = obj.validAxes();
            switch(opts.YTickLabels)
                case 'left'
                    ix = arrayfun(@(c)~c.UserData{1}.leftCol,validAxes);
                    set(validAxes(ix),'YTickLabel',[]);
            end
            switch(opts.XTickLabels)
                case 'bottom'
                    ix = arrayfun(@(c)~c.UserData{1}.botRow,obj.Axes(ishandle(obj.Axes)));
                    set(validAxes(ix),'xTickLabel',[]);
            end
            if size(obj.Axes,1)==size(obj.Axes,2)
                switch(opts.diagonal)
                    case 'remove'
                        delete(obj.Axes(logical(eye(size(obj.Axes)))));
                    case 'hide'
                        set(obj.Axes(logical(eye(size(obj.Axes)))),'Visible','off');
                end
            end
        end

        function [] = formatByLoc(obj,ax,opts)
            arguments
                obj
                ax = []
                opts.xTickLabels = 'bottom'
                opts.yTickLabels = 'left'
            end

            if iscell(ax)
                for i=1:numel(ax)
                    obj.formatByLoc(ax{i});
                end
                return
            end

            if isempty(ax)
                switch(opts.xTickLabels)
                    case 'bottom'
                        set(setdiff(obj.Axes(:),obj.getAxesByLoc("bottom")),'XTickLabel',[]);
                    case 'top'
                        set(setdiff(obj.Axes(:),obj.getAxesByLoc("top")),'XTickLabel',[]);
                end
                switch(opts.yTickLabels)
                    case 'left'
                        set(setdiff(obj.Axes(:),obj.getAxesByLoc("left")),'YTickLabel',[]);
                    case 'right'
                        set(setdiff(obj.Axes(:),obj.getAxesByLoc("right")),'YTickLabel',[]);
                end
            else
                locs = obj.getLocations(ax);
                switch(opts.xTickLabels)
                    case 'bottom'
                        set(setdiff(ax,ax(locs(:,2)=="bottom")),'XTickLabel',[]);
                    case 'top'
                        set(setdiff(ax,ax(locs(:,2)=="top")),'XTickLabel',[]);
                    case 'none'
                        set(ax,'XTickLabel',[]);
                end
                switch(opts.yTickLabels)
                    case 'left'
                        set(setdiff(ax,ax(locs(:,1)=="left")),'YTickLabel',[]);
                    case 'right'
                        set(setdiff(ax,ax(locs(:,1)=="right")),'YTickLabel',[]);
                    case 'none'                   
                        set(ax,'YTickLabel',[]);
                end
            end
        end

        function [locations] = getLocations(obj,ax)
            arguments
                obj
                ax = []
            end
            if isempty(ax), ax = obj.Axes; end
            axpos = vertcat(ax.Position);
            axxy = axpos(:,1:2)+axpos(:,3:4);
            corners = [0 0; 1 0; 1 1; 0 1]; %bl, br, tr, tl
            d = pdist2(corners,axxy)';
            [~,ixs] = min(d,[],1);
            locations = repmat("",numel(ax),2);
            locations(ixs(1),:) = ["left" "bottom"];
            locations(ixs(2),:) = ["right" "bottom"];
            locations(ixs(3),:) = ["right" "top"];
            locations(ixs(4),:) = ["left" "top"];
            tol = @(a,b)abs(a-b)<1e-3;
            locations(tol(axxy(:,1),axxy(ixs(1),1)),1) = "left";
            locations(tol(axxy(:,1),axxy(ixs(2),1)),1) = "right";
            locations(tol(axxy(:,2),axxy(ixs(1),2)),2) = "bottom";
            locations(tol(axxy(:,2),axxy(ixs(3),2)),2) = "top";            
        end

        function [ax] = getAxesByLoc(obj,varargin)
            axLoc = obj.getLocations();
            ixs = any(ismember(axLoc,string(varargin)),2);
            ax = obj.Axes(ixs);
        end

        function [] = identifyAxes(obj)
            for i=1:numel(obj.Axes)
                th(i) = text(mean(obj.Axes(i).XLim),mean(obj.Axes(i).YLim),num2str(i),...
                    'FontSize',36,'Parent',obj.Axes(i));
            end
            t = timer('StartDelay',10,'TimerFcn',@(~,~)delete(th),'ExecutionMode','singleShot');
            start(t);
        end

        function [] = labelAxes(obj,axObjs,labels,positions,opts)
            arguments
                obj
                axObjs
                labels
                positions = "nwo"
                opts.FontSize = 24
                opts.FontWeight = 'bold';
                opts.FontAngle = 'normal';
                opts.Interpreter = 'tex';
                opts.HorizontalAlignment = [];
                opts.VerticalAlignment = [];
                opts.prefix = '';
                opts.brackets = '';
                opts.Rotation = 0;
                opts.Offsets = [nan nan]
                opts.usePlotBox = true
            end
            if isempty(obj.BackgroundAxes)
                obj.createBackgroundAxes();
                obj.BackgroundAxes.Visible = 'off';
                uistack(obj.BackgroundAxes,'top');
            end
            positions = string(positions);

            if isempty(axObjs)
                axObjs = obj.Axes;
                switch(obj.AxesHandleForm)
                    case 'grid'
                        axObjs = axObjs';
                end
                axObjs = axObjs(ishandle(axObjs));
            end

            if ischar(labels) && any(strcmp(labels,{'letters' 'capitals' 'numbers' 'romannumerals'}))
                opts.prefix = labels;
                labels = "";
            end
            if isscalar(positions), positions = repmat(positions,size(axObjs)); end
            if isscalar(labels), labels = repmat(labels,size(axObjs)); end

            if ~isempty(opts.prefix)
                switch(opts.prefix)
                    case 'letters'
                        prefixes = arrayfun(@(c)string(char(c)),96+(1:numel(axObjs))');
                    case 'capitals'
                        prefixes = arrayfun(@(c)string(char(c)),64+(1:numel(axObjs))');
                    case 'numbers'
                        prefixes = string(1:numel(axObjs));
                    case 'romannumerals'
                        prefixes = arrayfun(@(c)string(toRoman(c)),(1:numel(axObjs))');

                end
                switch(opts.brackets)
                    case '()'
                        prefixes = "(" + prefixes + ")";
                    case ')'
                        prefixes = prefixes + ")";
                    case '[]'
                        prefixes = "[" + prefixes + "]";
                end
                labels = prefixes + labels(:);
            end

            labIxs = length(obj.AxesLabels)+(1:length(labels));
            for i=1:length(labIxs)
                alignInfo = stTextAlign.align(axObjs(i),positions(i),usePlotBox=opts.usePlotBox);
                obj.AxesLabels(labIxs(i)) = text(nan,nan,labels{i}, ...
                    FontSize=opts.FontSize,Parent=obj.BackgroundAxes,...
                    Interpreter=opts.Interpreter,Rotation=opts.Rotation,...
                    FontWeight=opts.FontWeight,FontAngle=opts.FontAngle);
                if ~isempty(opts.HorizontalAlignment)
                    alignInfo.HorizontalAlignment = opts.HorizontalAlignment;
                end
                if ~isempty(opts.VerticalAlignment)
                    alignInfo.VerticalAlignment = opts.VerticalAlignment;
                end
                stTextAlign.setAlignProps(obj.AxesLabels(labIxs(i)),alignInfo);
            end

            if ~all(isnan(opts.Offsets))
                shiftposx(obj.AxesLabels(labIxs),opts.Offsets(1));
                shiftposy(obj.AxesLabels(labIxs),opts.Offsets(2));
            end

        end

        function [h] = labelFigure(obj,str,pos,opts)
            arguments
                obj
                str
                pos
                opts.FontSize = 24
                opts.FontWeight = 'bold';
                opts.Interpreter = 'tex'
            end
            if isempty(obj.BackgroundAxes)
                obj.createBackgroundAxes();
                obj.BackgroundAxes.Visible = 'off';
            end
            alignInfo = stTextAlign.align(obj.BackgroundAxes,pos);
            h = text(alignInfo.Position(1),alignInfo.Position(2),str,...
                Parent=obj.BackgroundAxes,HorizontalAlignment=alignInfo.HorizontalAlignment,...
                VerticalAlignment=alignInfo.VerticalAlignment,FontSize=opts.FontSize,FontWeight=opts.FontWeight,...
                Interpreter=opts.Interpreter);
        end

        function [] = subAxes(obj,subAx,panelSpec,opts)

            arguments
                obj
                subAx
                panelSpec
                opts.intMargins = [0 0];
                opts.extMargins = [0 0 0 0];
                opts.Tags = []
            end

            panelLayout = obj.processPanelSpec(panelSpec);
            layoutProps = obj.layoutProperties(panelLayout);

            %calculate new positions in normalized coords
            newPositions = obj.panelLayoutToCoords(panelLayout,opts.extMargins,opts.intMargins);

            if isgraphics(subAx,'Axes')
                subAx = find(obj.Axes==subAx);
                oldTag = obj.Tags(subAx);
            else            
                oldTag = obj.Tags(subAx);
            end

            %transform to sub-axis positions
            subAxPositions = obj.Axes(subAx).Position;
            newPositions(:,1) = subAxPositions(1) + newPositions(:,1).*subAxPositions(3);
            newPositions(:,2) = subAxPositions(2) + newPositions(:,2).*subAxPositions(4);
            newPositions(:,3:4) = newPositions(:,3:4).*subAxPositions(3:4);

            %create new axes and remove substituted axis
            newAxes = gobjects(size(newPositions,1),1);
            newTags = repmat("",size(newPositions,1),1);
            for i=1:numel(newAxes)
                if i==1
                    newAxes(i) = obj.Axes(subAx);
                    set(newAxes(i),'Position',newPositions(i,:));
                else
                    newAxes(i) = axes(obj.Axes(subAx).Parent,'Position',newPositions(i,:), ...
                        'ActivePositionProperty','position');
                end
                if ~isempty(newAxes(i).Toolbar), newAxes(i).Toolbar.Visible = 'off'; end
                newAxes(i).UserData{1} = layoutProps(i);
                newAxes(i).UserData{2} = obj;  
                if ~isempty(opts.Tags)
                    newTags(i) = opts.Tags(i);
                else
                    newTags(i) = oldTag;
                end
            end
            % Retain the original target axes as the first panel.
            switch(obj.AxesHandleForm)
                case 'tensor'
                    if ndims(obj.Axes)~=3
                        ax = repmat(gobjects(obj.PanelSpec),[1 1 numel(newAxes)]);
                        for a=1:size(obj.Axes,1)
                            for b=1:size(obj.Axes,2)
                                ax(a,b,1) = obj.Axes(a,b,1);
                            end
                        end
                        obj.Axes = ax;
                    end
                    [r,c] = ind2sub(obj.PanelSpec(1:2),subAx);
                    obj.Axes(r,c,:) = newAxes;
                otherwise
                    obj.Axes = [obj.Axes(1:subAx-1); newAxes; obj.Axes(subAx+1:end)]; % Update the Axes property with new axes
                    obj.Tags = [obj.Tags(1:subAx-1); newTags; obj.Tags(subAx+1:end)];
            end

            obj.currentAxes = newAxes(1);
            obj.CurrentAxesIx = subAx;
            set(obj.Handle,'CurrentAxes',obj.currentAxes);
        end

        function [] = arraySubAxes(obj,axesObjs,panelLayout,opts)
            arguments
                obj
                axesObjs
                panelLayout
                opts.extMarg = [0 0 0 0]
                opts.intMarg = [0 0]
                opts.Tags = []
            end
            currAxesIx = obj.CurrentAxesIx;
            for i=1:numel(axesObjs)
                obj.subAxes(axesObjs(i),panelLayout,extMargins=opts.extMarg,intMargins=opts.intMarg,Tags=opts.Tags);                
            end
            obj.CurrentAxesIx = currAxesIx;
            obj.currentAxes = obj.Axes(obj.CurrentAxesIx);
            axes(obj.currentAxes);
        end

        function [layoutProps] = layoutProperties(~,panelLayout)
            leftCol = panelLayout==panelLayout(:,1) & ~isnan(panelLayout);
            rightCol = panelLayout==panelLayout(:,end) & ~isnan(panelLayout);
            topRow = panelLayout==panelLayout(1,:) & ~isnan(panelLayout);
            botRow = panelLayout==panelLayout(end,:) & ~isnan(panelLayout);
            panNums = unique(panelLayout(~isnan(panelLayout)));
            for i=panNums(:)'
                layoutProps(i).leftCol = leftCol(find(panelLayout==i,1,'first'));
                layoutProps(i).rightCol = rightCol(find(panelLayout==i,1,'first'));
                layoutProps(i).topRow = topRow(find(panelLayout==i,1,'first'));
                layoutProps(i).botRow = botRow(find(panelLayout==i,1,'first'));
            end
        end

        function [] = setFontSize(obj,opts)
            arguments
                obj
                opts.axes = []
                opts.labels = []
                opts.xlabel = []
                opts.ylabel = []
                opts.axis = []
                opts.title = []
                opts.legend = []
                opts.text = []
            end

            if isempty(opts.axes)
                opts.axes = obj.Axes;
            end
            ax = [opts.axes(:); obj.BackgroundAxes];
            ax = ax(ishandle(ax));

            CH = arrayfun(@(c){vertcat(allchild(c))},ax);
            CH = [ax; vertcat(CH{:})];

            O = table(CH,'VariableNames',{'object'});
            O.type = arrayfun(@(c){get(c,'Type')},O.object);

            ax_ch = {'axis','labels','ylabel','xlabel','text','title','legend'};
            for i=1:height(O)
                switch(O.type{i})
                    case 'axes'
                        for j=1:length(ax_ch)
                            if ~isempty(opts.(ax_ch{j}))
                                switch(ax_ch{j})
                                    case {'labels','label'}
                                        O.object(i).XLabel.FontSize = opts.labels;
                                        O.object(i).YLabel.FontSize = opts.labels;
                                        O.object(i).ZLabel.FontSize = opts.labels;
                                    case 'xlabel'
                                        O.object(i).XLabel.FontSize = opts.xlabel;
                                    case 'ylabel'
                                        O.object(i).YLabel.FontSize = opts.ylabel;
                                    case 'title'
                                        O.object(i).Title.FontSize = opts.title;
                                    case {'axis','axes'}
                                        xlim(O.object(i),O.object(i).XLim);
                                        ylim(O.object(i),O.object(i).YLim);
                                        O.object(i).FontSize = opts.axis;
                                    case 'legend'
                                        if ~isempty(O.object(i).Legend)
                                            O.object(i).Legend.FontSize = opts.legend;
                                        end
                                end
                            end
                        end
                    case 'text'
                        if ~isempty(opts.text)
                            O.object(i).FontSize = opts.text;
                        end
                end

            end
        end

        function [] = centerAxes(obj,ax,dim)
            arguments
                obj
                ax = []
                dim = 'both'
            end
            if isempty(ax), ax = obj.Axes; end
            for i=1:numel(ax)
                obj.centerAxis(ax(i),dim);
            end
        end

        function [] = centerAxis(obj,ax,dim)
            dataLims = obj.getVisibleObjectLims(ax);
            if ~isequal(size(dataLims),[2 2])
                return  % nothing visible to center against
            end
            axisLims = [ax.XLim; ax.YLim];
            shifts = mean(dataLims,2) - mean(axisLims,2);
            switch(lower(dim))
                case {'horizontal','x'}
                    ax.XLim = ax.XLim + shifts(1);
                case {'vertical','y'}
                    ax.YLim = ax.YLim + shifts(2);
                case {'both','xy'}
                    ax.XLim = ax.XLim + shifts(1);
                    ax.YLim = ax.YLim + shifts(2);
            end
        end

        function [lims] = getVisibleObjectLims(~,ax)
            ch = get(ax,'Children');
            vis = string(get(ch,'Visible'));
            ch = ch(vis=="on");
            xd = []; yd=[];
            for i=1:length(ch)
                switch(ch(i).Type)
                    case 'text'
                        xd = [xd; ch(i).Extent(1); sum(ch(i).Extent([1 3]))];
                        yd = [yd; ch(i).Extent(2); sum(ch(i).Extent([2 4]))];
                    otherwise
                        xd = [xd; ch(i).XData(:)];
                        yd = [yd; ch(i).YData(:)];
                end
            end

            lims = [minmax(xd'); minmax(yd')];
        end

        function [] = copyToClipboard(obj)
            obj.copy();
        end

        function [] = setTickLen(obj,tickLen)
            arguments
                obj
                tickLen = 0.0025
            end
            set(obj.Axes(:),'TickLen',tickLen*[1 1]);
        end

        function [] = rescaleAxes(obj,ax,dim,scalevals1,scalevals2)
            arguments
                obj
                ax = []
                dim = 'both'
                scalevals1 = [0 0]
                scalevals2 = [0 0]
            end
            if isempty(ax), ax = obj.Axes; end
            if ~ishandle(ax), ax = obj.Axes(ax); end
            ax = ax(:);
            for i=1:length(ax)
                switch(dim)
                    case 'x'
                        stFig.rescaleaxes(ax(i),'x',scalevals1);
                    case 'y'
                        stFig.rescaleaxes(ax(i),'y',scalevals1);
                    case {'both','xy'}
                        stFig.rescaleaxes(ax(i),'x',scalevals1);
                        stFig.rescaleaxes(ax(i),'y',scalevals2);
                end
            end
        end
        function [] = equateGridSpacing(obj,ax,opts)
            arguments
                obj
                ax = []
                opts.spacing = []
            end
            if isempty(ax), ax = obj.Axes; end
            for i=1:numel(ax)
                spacing = opts.spacing;
                xticks = ax(i).XTick;
                yticks = ax(i).YTick;
                xts = mode(diff(xticks));
                yts = mode(diff(yticks));
                if isempty(spacing)
                    spacing = max([xts yts]);
                end
                ax(i).XTick = xticks(1):spacing:xticks(end);
                ax(i).YTick = yticks(1):spacing:yticks(end);
            end
        end

        function [] = Ax(obj,num)            
            if isscalar(num)
                if isnumeric(num)
                    axes(obj.Axes(num));
                    obj.CurrentAxesIx = num;
                    obj.currentAxes = obj.Axes(num);
                elseif ishandle(num)
                    axes(num);                    
                    obj.CurrentAxesIx = find(obj.Axes==num);
                    obj.currentAxes = num;
                end
            elseif numel(num)==2
                axes(obj.Axes(num(1),num(2)));
            elseif numel(num)==3
                axes(obj.Axes(num(1),num(2),num(3)));
            end
        end

        function [ax] =  getTagged(obj,tags,ixs)            
            ax = obj.Axes(ismember(obj.Tags,tags));
            if nargin==3
                ax = ax(ixs);
            end
        end

        function [] = next(obj)
            if isempty(obj.CurrentAxesIx)
                obj.CurrentAxesIx = 1;
                set(obj.Handle,'CurrentAxes',obj.Axes(obj.CurrentAxesIx));
                obj.currentAxes = obj.Axes(obj.CurrentAxesIx);                
                return
            end
            
            Naxes = numel(obj.Axes(ishandle(obj.Axes)));
            nextIx = obj.CurrentAxesIx+1;
            if nextIx>max(obj.PanelLayout(:))
                nextIx = 1;
            end
            [r,c] = find(obj.PanelLayout==nextIx);
            
            obj.CurrentAxesIx =min(obj.CurrentAxesIx+1,Naxes);
            switch(obj.AxesHandleForm)
                case 'grid'
                    set(obj.Handle,'CurrentAxes',obj.Axes(r,c));
                    obj.currentAxes = obj.Axes(r,c);
                otherwise
                    set(obj.Handle,'CurrentAxes',obj.Axes(obj.CurrentAxesIx));
                    obj.currentAxes = obj.Axes(obj.CurrentAxesIx);
            end
            
        end

        function [] = optimizeFontSize(obj,ch,opts)
            arguments
                obj
                ch = []
                opts.constraints = ["axes" "figure" "data" "pixels"]
                opts.DataProportion = 0.05
                opts.FigProportion = 0.03
                opts.MaxCharWidthPx = 14
                opts.Dimension = 'y'
                opts.Scale = 1
            end
            % ---- Resolve figure and axes list
            if ishghandle(h,'figure')
                fig = h;
                axs = findall(fig,'Type','axes');
            elseif ishghandle(h,'axes')
                axs = h;
                fig = ancestor(h,'figure');
            else
                error('h must be a figure or axes handle.');
            end
            axs = axs(isvalid(axs) & strcmp(get(axs,'Visible'),'on')); %#ok<GSET>

            if isempty(axs)
                return;
            end

            % ---- Figure height in pixels
            figUnits = fig.Units;
            fig.Units = 'pixels';
            figPos   = fig.Position;
            figHpx   = max(1, figPos(4));
            fig.Units = figUnits;

            % ---- Pixels-per-point for conversion
            ppp = get(groot,'ScreenPixelsPerInch')/72;  % pixels per typographic point

            % ---- Collect all text-like objects per axes
            for ax = reshape(axs,1,[])
                if ~isvalid(ax), continue; end

                % Plotbox size (pixels)
                axUnits = ax.Units;
                ax.Units = 'pixels';
                if isprop(ax,'InnerPosition')
                    pb = ax.InnerPosition;      % [x y w h] in px
                else
                    pb = ax.Position;           % fallback
                end
                ax.Units = axUnits;

                plotW = max(1, pb(3));
                plotH = max(1, pb(4));
                switch dim
                    case 'y',   dataPx = opt.DataProportion * plotH;
                    case 'x',   dataPx = opt.DataProportion * plotW;
                    case 'min', dataPx = opt.DataProportion * min(plotW, plotH);
                end

                figPx = opt.FigProportion * figHpx;

                % Build list: plain text + title/xlabel/ylabel/zlabel
                txt = findall(ax, 'Type','text');

                deco = struct('obj',[]); idx=0;
                % Titles/labels are distinct classes; we include if present
                if isprop(ax,'Title')   && ~isempty(ax.Title);   idx=idx+1; deco(idx).obj=ax.Title;   end
                if isprop(ax,'XLabel')  && ~isempty(ax.XLabel);  idx=idx+1; deco(idx).obj=ax.XLabel;  end
                if isprop(ax,'YLabel')  && ~isempty(ax.YLabel);  idx=idx+1; deco(idx).obj=ax.YLabel;  end
                if isprop(ax,'ZLabel')  && ~isempty(ax.ZLabel);  idx=idx+1; deco(idx).obj=ax.ZLabel;  end

                % Optionally ticks (rulers)
                rulers = {};
                if opt.IncludeTicks
                    if isprop(ax,'XAxis'), rulers{end+1} = ax.XAxis; end
                    if isprop(ax,'YAxis'), rulers{end+1} = ax.YAxis; end
                    if isprop(ax,'ZAxis'), rulers{end+1} = ax.ZAxis; end
                end

                % Apply to raw text objects
                for t = reshape(txt,1,[])
                    if ~isvalid(t), continue; end
                    sPx = computeCapForObject(t, figPx, dataPx, opt.MaxCharWidthPx, ax);
                    setFontSizePixels(t, sPx * opt.Scale, ppp);
                end

                % Apply to decorations (Title/Labels)
                for k = 1:numel(deco)
                    o = deco(k).obj;
                    if ~isvalid(o), continue; end
                    sPx = computeCapForObject(o, figPx, dataPx, opt.MaxCharWidthPx, ax);
                    setFontSizePixels(o, sPx * opt.Scale, ppp);
                end

                % Apply to ticks if requested (rulers don't have Interpreter, but have FontName/Weight)
                for r = rulers
                    rr = r{1};
                    if ~isvalid(rr), continue; end
                    % Fake a text-like proxy using ruler font props to estimate char width
                    proxy.FontName   = getProp(rr,'FontName','');
                    proxy.FontWeight = getProp(rr,'FontWeight','normal');
                    proxy.FontAngle  = getProp(rr,'FontAngle','normal');
                    proxy.Interpreter = 'tex'; % ticks use plain rendering
                    sPx = min([dataPx, figPx, charWidthCapPx(proxy, opt.MaxCharWidthPx, ax)]);
                    % Rulers expect points directly
                    set(rr, 'FontSize', (sPx * opt.Scale) / ppp);
                end
            end
        end

        function [] = moveWindow(obj,screenshift,opts)
            arguments
                obj
                screenshift = 1
                opts.maximize = true
            end
            shiftposx(obj.Handle,screenshift);
            if opts.maximize, set(obj.Handle,'WindowState','maximized'); end
        end

        function [ax] = validAxes(obj,axin)
            if nargin==2 && ~isempty(axin)
                ax = axin(ishandle(axin));
            else
                ax = obj.Axes(ishandle(obj.Axes));
            end
        end

        function [] = Tight(obj,ax)
            arguments
                obj
                ax = []
            end
            if isempty(ax), ax=obj.Axes; end
            obj.tight(ax);
        end
        function [] = EquateAxesSizes(obj,ax,dims)
            arguments
                obj
                ax = []
                dims = 'xy'
            end            
            obj.equateAxesSizes(obj.validAxes(ax),dims=dims);
        end
        function [] = SymmetrizeLims(obj,ax,dims)
            arguments
                obj
                ax = []
                dims = 'xy'
            end            
            obj.symmetrizeLims(obj.validAxes(ax),dims);
        end
        function [] = SetAxesProps(obj,ax,props)
            arguments
                obj
                ax = []
                props = {}
            end
            obj.setAxesProps(obj.validAxes(ax),props);
        end
        function [] = RescaleAxes(obj,ax,dims,varargin)
            arguments
                obj
                ax = []
                dims='xy'
            end
            arguments (Repeating)
                varargin
            end
            obj.rescaleAxes(obj.validAxes(ax),dims,varargin{:});
        end
        function [lh] = objConnectionLine(obj,objs,pos,opts)
            arguments
                obj
                objs
                pos
                opts.color = [0 0 0]
                opts.lineWidth = 2
                opts.lineStyle = '-'
            end

            for i=1:numel(objs)
                switch(objs(i).Type)
                    case 'axes'
                        posxy(i,:) = objpos(objs(i),pos(i));
                        
                    otherwise
                        posxy(i,:) = obj.data2nfu(objs(i).Parent,objpos(objs(i),pos(i)));

                end
            end
            
            if isempty(obj.BackgroundAxes)
                obj.createBackgroundAxes();
            end
            lh = plot(posxy(:,1),posxy(:,2),Color=opts.color,LineStyle=opts.lineStyle,LineWidth=opts.lineWidth,...
                Parent=obj.BackgroundAxes);

        end

        function [] = assignTags(obj)
            tags = unique(obj.Tags);
            for i=1:numel(tags)
                ax = obj.getTagged(tags(i));                
                assignin("caller","ax"+tags(i),ax);
            end
        end

    end

    methods (Static)

        function tf = isLiveScript()
            % Best-effort detection of whether the CALLING code is
            % executing inside a Live Script run via the Live Editor.
            % There is no documented, stable public API for this as of
            % R2026a. Live Editor execution never puts the .mlx path in
            % the call stack — instead the outermost frame is a temp file
            % named LiveEditorEvaluationHelper*.m under a Temp\Editor_*
            % folder (confirmed empirically via breakpoint on R2026a
            % Windows); that's the signal checked here. dbstack is a
            % standard, documented function, so only the inference from
            % its contents is a heuristic, and this always falls back to
            % false — i.e. desktop/plain-script behavior — rather than
            % erroring if that's inconclusive.
            %
            % If auto-detection guesses wrong for your setup, override it
            % explicitly for the CURRENT PROCESS ONLY. Deliberately an
            % environment variable rather than setpref: setpref persists
            % to disk and is shared by every MATLAB process for this
            % user, so a crashed/interrupted run that sets it true and
            % never gets to undo it silently makes every figure --
            % including in the user's own later interactive session --
            % behave as if it's in a Live Script (this happened during
            % development: an interrupted export left LiveScriptOverride
            % stuck on, and every desktop-mode figure after that silently
            % went offscreen instead of appearing normally). An env var
            % only exists for the process it's set on, so there's nothing
            % to leave stuck:
            %   setenv('STFIG_LIVESCRIPT_OVERRIDE','1')  % force live-script behavior
            %   setenv('STFIG_LIVESCRIPT_OVERRIDE','0')  % force desktop behavior
            %   setenv('STFIG_LIVESCRIPT_OVERRIDE','')   % return to auto-detection
            ov = getenv('STFIG_LIVESCRIPT_OVERRIDE');
            if any(strcmp(ov,{'1','0'}))
                tf = strcmp(ov,'1');
                return;
            end
            tf = false;
            try
                stack = dbstack('-completenames');
                files = {stack.file};
                tf = any(contains(files,'LiveEditorEvaluationHelper','IgnoreCase',true)) || ...
                     any(endsWith(files,'.mlx','IgnoreCase',true));
            catch
                tf = false;
            end
        end

        function figPt = data2nfu(ax, xy)


            % Ensure 2D
            dataPt = xy(1:2);
            x = dataPt(1);
            y = dataPt(2);

            % Store and temporarily change units to normalized
            fig = ancestor(ax, 'figure');
            oldFigUnits = fig.Units;
            oldAxUnits  = ax.Units;
            fig.Units   = 'normalized';
            ax.Units    = 'normalized';

            % --- Get plot-box rectangle in normalized figure units ---
            % Prefer tightPosition (plotting area == plot box)
            try
                pb = tightPosition(ax);      % [left bottom width height]
            catch
                % Fallback: use axes Position (less accurate for axis equal/image)
                pb = ax.Position;
            end

            % Data limits
            xlim = ax.XLim;
            ylim = ax.YLim;

            % Normalized within data range
            xn = (x - xlim(1)) / (xlim(2) - xlim(1));
            yn = (y - ylim(1)) / (ylim(2) - ylim(1));

            % Map into figure normalized using the *plot box* rectangle
            figX = pb(1) + pb(3) * xn;
            figY = pb(2) + pb(4) * yn;

            figPt = [figX figY];

            % Restore units
            fig.Units = oldFigUnits;
            ax.Units  = oldAxUnits;


        end


        function [lh] = drawSpanLine(ax,opts)
            arguments
                ax
                opts.dims = 'y';
                opts.color = [0 0 0];
                opts.lineStyle = ':';
                opts.lineWidth = 1;
                opts.stack = 'bottom'
                opts.value = 0;
            end

            props = {'lineStyle',opts.lineStyle,'lineWidth',opts.lineWidth,'color',opts.color};

            values = opts.value*[1 1];

            for i=1:numel(ax)
                switch(opts.dims)
                    case 'x'
                        lh = plot(ax(i).XLim,values,'-',props{:},'parent',ax(i));
                    case 'y'
                        lh = plot(values,ax(i).YLim,'-',props{:},'parent',ax(i));
                    case 'xy'
                        lh(1) = plot(ax(i).XLim,values,'-',props{:},'parent',ax(i));
                        lh(2) = plot(values,ax(i).YLim,'-',props{:},'parent',ax(i));
                end
                try
                    uistack(lh,opts.stack);
                catch
                    %error for multi-axis plots
                end
            end

        end

        function [cmap] = colormap(ncolors,varargin)
            arguments
                ncolors = 256
            end
            arguments (Repeating)
                varargin
            end

            c = varargin;
            cmap = [];

            ls = @(col1,col2,dim)linspace(col1(dim),col2(dim),ncolors)';
            for i=1:length(varargin)-1
                cmap = [cmap; ls(c{i},c{i+1},1) ls(c{i},c{i+1},2) ls(c{i},c{i+1},3)];
            end
            ix = find(all(diff(cmap)==0,2));
            cmap = cmap(setdiff(1:size(cmap,1),ix),:);
        end

        function [panelLayout,panelsStr] = processPanelSpec(panelSpec)
            if isscalar(panelSpec) %make as square as possible
                nr = floor(sqrt(panelSpec));
                nc = ceil(panelSpec/nr);
                panelLayout = reshape(1:nr*nc,nc,[])';
                panelLayout(panelLayout>panelSpec) = nan;
            elseif all(size(panelSpec)==[1 2])
                panelLayout = reshape(1:prod(panelSpec),panelSpec(2),[])';
            else
                panelLayout = panelSpec;
            end
            panelsStr = num2str(panelLayout);
        end

        function [panelPositions] = panelLayoutToCoords(panelLayout,exteriorMargins,interiorMargins)
            xmarg = exteriorMargins([1 3]);
            ymarg = exteriorMargins([2 4]);
            imarg = interiorMargins;

            axn = flipud(panelLayout);
            nr = size(axn,1);
            nc = size(axn,2);
            arng = 1 - sum(xmarg) + imarg(1);
            brng = 1 - sum(ymarg) + imarg(2);

            %convert to offsets
            w = arng/nc;
            h = brng/nr;

            xoff = xmarg(1):w:(1-xmarg(2));
            yoff = ymarg(1):h:(1-ymarg(2));

            if length(xoff)<size(axn,2)
                xoff = linspace(xmarg(1),(1-xmarg(2)),size(axn,2)); w = mean(diff(xoff));
            end
            if length(yoff)<size(axn,1)
                yoff = linspace(ymarg(1),(1-ymarg(2)),size(axn,1)); h = mean(diff(yoff));
            end

            nax = length(unique(axn(~isnan(axn(:)))));
            axposc = cell(1,nax);
            for y=1:nr
                for x=1:nc
                    if isnan(axn(y,x)), continue; end
                    axposc{axn(y,x)} = [axposc{axn(y,x)}; xoff(x) yoff(y) xoff(x)+(w-imarg(1)) yoff(y)+(h-imarg(2))];
                end
            end

            for j=1:length(axposc)
                panelPositions(j,:) = [min(axposc{j}(:,1)) min(axposc{j}(:,2)) max(axposc{j}(:,3)) max(axposc{j}(:,4))]; %#ok<AGROW>

            end

            panelPositions(:,3:4) = panelPositions(:,3:4)-panelPositions(:,1:2);
        end

        function rescaleaxes(ax,dim,ext,ext2)
            ax = ax(:);
            if isscalar(ext), ext = repmat(ext,1,2); end
            for i=1:numel(ax)
                pos = get(ax(i),'Position');
                hold(ax(i),'on');
                switch(dim)
                    case 'x'
                        ax(i).XLim = ax(i).XLim + [-1 1]*diff(ax(i).XLim).*ext;
                    case 'y'
                        ax(i).YLim = ax(i).YLim + [-1 1]*diff(ax(i).YLim).*ext;
                    case 'xy'
                        if isscalar(ext2), ext2 = repmat(ext2,1,2); end
                        ax(i).XLim = ax(i).XLim + [-1 1]*diff(ax(i).XLim).*ext;
                        ax(i).YLim = ax(i).YLim + [-1 1]*diff(ax(i).YLim).*ext2;
                end
                set(ax(i),'Position',pos);
            end
        end

        function equateAxesSizes(ax,opts)
            arguments
                ax
                opts.dims = 'xy'
            end
            ax = ax(:);
            ax = ax(ishandle(ax));
            if numel(ax)<=1, return; end

            if contains(opts.dims,'x')
                xlims = cell2mat(get(ax,'Xlim'));
                xlims = minmax(xlims(:)');
                set(ax,'XLim',xlims);
            end
            if contains(opts.dims,'y')
                ylims = cell2mat(get(ax,'Ylim'));
                ylims = minmax(ylims(:)');
                set(ax,'YLim',ylims);
            end
        end
        function equateAxesScales(ax)
            arguments
                ax
            end
            ax = ax(ishandle(ax));
            if numel(ax)<=1, return; end

            xlims = cell2mat(get(ax,'Xlim'));
            ylims = cell2mat(get(ax,'Ylim'));
            dx = diff(xlims,[],2);
            dy = diff(ylims,[],2);
            scx = abs(dx);
            scy = abs(dy);

            ix = scx<scy;
            xlims(ix,:) = mean(xlims(ix,:),2) + scy(ix)*[-1 1]/2;
            ix = scx>scy;
            ylims(ix,:) = mean(ylims(ix,:),2) + scx(ix)*[-1 1]/2;

            for i=1:numel(ax)
                set(ax(i),'xlim',xlims(i,:),'ylim',ylims(i,:));
            end

        end

        function symmetrizeLims(ax,dims)
            arguments
                ax = []
                dims = 'y'
            end
            ax = ax(ishandle(ax));
            for i=1:numel(ax)
                switch(dims)
                    case 'x'
                        set(ax(i),'xlim',max(abs(xlim(ax(i))))*[-1 1]);
                    case 'y'
                        set(ax(i),'ylim',max(abs(ylim(ax(i))))*[-1 1]);
                    case 'xy'
                        set(ax(i),'xlim',max(abs(xlim(ax(i))))*[-1 1]);
                        set(ax(i),'ylim',max(abs(ylim(ax(i))))*[-1 1]);
                end
            end
        end

        function tight(ax)
            axis(ax(ishandle(ax)),'tight');
        end

        function setAxesProps(ax,props)
            set(ax(ishandle(ax)),props{:});
        end

        function shiftAxes(ax,dim,shift)
            ax = ax(ishandle(ax));
            for i=1:numel(ax)
                switch(dim)
                    case 'x'
                        ax(i).Position(1) = ax(i).Position(1)+shift;
                    case 'y'
                        ax(i).Position(2) = ax(i).Position(2)+shift;
                end
            end
        end

        function toggleVisible(objs)
            for i=1:length(objs)
                switch(ishandle(objs(i)))
                    case true
                        switch(get(objs(i),'Visible'))
                            case 'on'
                                set(objs(i),'Visible','off');
                            case 'off'
                                set(objs(i),'Visible','on');
                        end
                end
            end
        end

        function layout = panelLayout(n,opts)
            arguments
                n
                opts.nrows = []
                opts.ncols = []
                opts.fillDir = 'horizontal'
            end
            if isempty(opts.nrows) & isempty(opts.ncols)
                nr = ceil(sqrt(n));
                nc = nr;
            elseif ~isempty(opts.nrows)
                nr = opts.nrows;
                nc = ceil(n/nr);
            else
                nc = opts.ncols;
                nr = ceil(n/nc);
            end
            switch(opts.fillDir)
                case 'horizontal'
                    layout = reshape(1:nr*nc,[],nr);
                    layout((n+1):end) = nan;
                    layout = layout';
                case 'vertical'
                    layout = reshape(1:nr*nc,[],nc);
                    layout((n+1):end) = nan;
            end

        end

        function [tile] = defTile(axpan,extMargins,intMargins,tags)
            arguments
                axpan = 1
                extMargins = [0 0]
                intMargins = [0 0]
                tags = ""
            end
            tile.panSpec = axpan;
            tile.extMargins = extMargins;
            tile.intMargins = intMargins;
            tile.tags = tags;
        end

        function [obj] = tileLayout(tile,varargin)
            obj = stFig(varargin{:});
            
            i = 1;
            while i<=numel(obj.Axes)            
                obj.subAxes(i,tile.panSpec,extMargins=tile.extMargins,intMargins=tile.intMargins);
                subIxs = i+(0:max(tile.panSpec)-1);
                for j=1:numel(subIxs)
                    obj.Axes(subIxs(j)).Tag = tile.tags(j);
                end
                i = max(subIxs)+1;
            end

            tags = string({obj.Axes.Tag});
            utags = unique(tags,'stable');
            %{
            for i=1:numel(utags)
                obj.Tags.(utags(i)) = obj.Axes(tags==utags(i));
            end
            %}
            obj.Tags = tags;
        end

        function [] = shiftTickValues(ax,dim,fcn,opts)
            arguments
                ax
                dim
                fcn
                opts.format = '%f'
            end

            for i=1:numel(ax)
                switch(dim)
                    case 'x'
                        xticks = ax(i).XTick;
                        xticks = fcn(xticks);
                        ax(i).XTickLabel = arrayfun(@(c){sprintf(opts.format,c)},xticks);

                    case 'y'
                        yticks = ax(i).YTick;
                        yticks = fcn(yticks);
                        ax(i).XTickLabel = arrayfun(@(c){sprintf(opts.format,c)},yticks);

                end
            end
        end

        function [] = moveYLabel(ax,factor)
            for i=1:numel(ax)
                ax(i).YLabel.Position(1) = ax(i).YLabel.Position(1) + factor*diff(xlim(ax(i)));
            end
        end
        function [] = moveXLabel(ax,factor)
            for i=1:numel(ax)
                ax(i).XLabel.Position(2) = ax(i).XLabel.Position(2) + factor*diff(ylim(ax(i)));
            end
        end
        function btn = wait(figh, opts)

            %       btn = 1 (left), 2 (middle), 3 (right)
            %       btn = nan if timeout expires without click.

            arguments
                figh = []
                opts.timeOut (1,1) double = inf
            end

            if isempty(figh) || ~ishandle(figh) || ~strcmp(get(figh,'Type'),'figure')
                figh = gcf;
            end

            % Shared variable for callback
            clickedBtn = nan;


            % Map SelectionType -> numeric button
            function b = mapSelectionType(sel)
                switch sel
                    case 'normal', b = 1;   % left
                    case 'extend', b = 2;   % middle (or shift-click depending on platform)
                    case 'alt',    b = 3;   % right
                    otherwise,     b = [];  % 'open' (double-click) or unknown
                end
            end

            % Mouse callback (ignore event.Button; use figure.SelectionType)
            function mouseDown(src, ~)
                clickedBtn = mapSelectionType(get(src,'SelectionType'));
                uiresume(figh);
            end
            % Attach temporary callback
            set(figh, 'WindowButtonDownFcn', @mouseDown);

            % Create timer for timeout (if finite)
            if isfinite(opts.timeOut) && opts.timeOut > 0
                t = timer('StartDelay', opts.timeOut, ...
                    'TimerFcn', @(~,~) uiresume(figh), ...
                    'ExecutionMode', 'singleShot');
                start(t);
            else
                t = [];
            end

            % Wait
            uiwait(figh);

            % Clean up
            if ~isempty(t) && isvalid(t)
                stop(t);
                delete(t);
            end
            set(figh, 'WindowButtonDownFcn', '');

            % Return result
            btn = clickedBtn;
        end

        function [varargout] = drawGridLines(ax,dims,opts)
            arguments
                ax
                dims = 'xy'
                opts.Color = [0.75*[1 1 1] 0.25]
                opts.LineWidth = 0.1
            end
            varargout = {};
            for i=1:numel(ax)
                if contains(dims,'x')
                    xt = ax(i).XTick;
                    gh = plot(repmat(xt,2,1),repmat(ax(i).YLim',1,numel(xt)),'-', ...
                        'color',opts.Color,'LineWidth',opts.LineWidth,Parent=ax(i));
                    uistack(gh,'top');
                    varargout{end+1} = gh;
                end
                if contains(dims,'y')
                    yt = ax(i).YTick;
                    gh = plot(repmat(ax(i).XLim',1,numel(yt)),repmat(yt,2,1),'-', ...
                        'color',opts.Color,'LineWidth',opts.LineWidth,Parent=ax(i));
                    uistack(gh,'top');
                    varargout{end+1} = gh;
                end
            end
            
        end

        function [h] = DrawAnnotation(fig,endpoint,label,alignment,opts)
            arguments
                fig
                endpoint
                label
                alignment                
                opts.lineType = 'wavy'
                opts.waveAmp = 0.1
                opts.fontName = 'calibri'
                opts.fontSize = 16
                opts.lineLength = 0.1
                opts.color = 0.25*[1 1 1];
                opts.Overlay = [] % newest overlay by default
            end
            
            previousAxes = fig.Handle.CurrentAxes;
            restore = onCleanup(@()set(fig.Handle,'CurrentAxes',previousAxes)); %#ok<NASGU>
            overlay = fig.getOverlay(opts.Overlay);
            %axes units
            if iscell(endpoint)
                endpoint = stFig.data2nfu(endpoint{1},endpoint{2});
                angle = alignspec2angle(alignment);
                startpoint = endpoint+[cos(angle) sin(angle)]*opts.lineLength;
            else %figure units
                angle = alignspec2angle(alignment);
                startpoint = endpoint+[cos(angle) sin(angle)]*opts.lineLength;
            end

            PP = [startpoint; endpoint];

            h.line = plot(overlay,PP(:,1),PP(:,2),'color',opts.color);
            switch(opts.lineType)
                case 'wavy'
                    draw_wavy_line(h.line,opts.waveAmp,nPoints=1000,freq=2);
            end
            h.arrow = draw_arrow(h.line,BaseLength=6,Anchor='base');
            align = stTextAlign.align(h.line,alignment+"o");
            h.text = text(overlay,h.line.XData(1),h.line.YData(1),label,...
                'HorizontalAlignment',align.HorizontalAlignment,...
                'VerticalAlignment',align.VerticalAlignment,...
                'FontSize',opts.fontSize,'FontName',opts.fontName);

        end
        function [] = hideTickLabel(ax,dim,num)
            arguments
                ax
                dim
                num = inf
            end
            for i=1:numel(ax)
                switch(dim)
                    case 'y'                        
                        ax(i).YTickLabelMode = 'auto'; 
                        ticklabels = ax(i).YTickLabel;
                        switch(num)
                            case inf
                                ticklabels{end} = '';
                            otherwise
                                ticklabels{1} = '';
                        end
                        ax(i).YTickLabel = ticklabels;
                    case 'x'                        
                        ax(i).XTickLabelMode = 'auto';
                        ticklabels = ax(i).XTickLabel;
                        switch(num)
                            case inf
                                ticklabels{end} = '';
                            otherwise
                                ticklabels{1} = '';
                        end
                        ax(i).XTickLabel = ticklabels;
                end

            end

        end
        function [] = matchAxisXY(ax,opts)
            arguments
                ax
                opts.all = false
            end

            xlims = cell2mat(get(ax,'XLim'));
            ylims = cell2mat(get(ax,'YLim'));
            lims = [xlims ylims];

            switch(opts.all)
                case false
                    lims = [min(lims,[],2) max(lims,[],2)];
                    for i=1:numel(ax)
                        set(ax(i),'XLim',lims(i,:),'YLim',lims(i,:));
                    end
                case true

            end
        end
        function [] = drawDiagonal(ax,opts)
            arguments
                ax
                opts.fcn = @(x)x;
                opts.color = 'k'
                opts.lineWidth = 1
                opts.lineStyle = '--';
            end
            for i=1:numel(ax)
                plot(ax(i).XLim,opts.fcn(ax(i).XLim), ...
                    LineStyle=opts.lineStyle, ...
                    Color=opts.color,LineWidth=opts.lineWidth,Parent=ax(i));
            end
        end
    end

end