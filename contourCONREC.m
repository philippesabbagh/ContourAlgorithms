% Philippe Sabbagh
% Draper
%
% CONREC subroutine by Paul Bourke July 1987 http://paulbourke.net/papers/conrec/

function [C] = contourPhilippeAdvanced(X,Y,Z,LEVELS)
plotFlag = 0;
plotTriangles = 0;
closeContours = 0; %EXPERIMENTAL
C = [];
[rows, cols] = size(X);
% plot(X,Y,'k.');
% Y = flip(Y);
figure; grid on; axis([min(X(:)) max(X(:)) min(Y(:)) max(Y(:))]); hold on;
%% SECTION 1: CONREC SUBROUTINE
for level = LEVELS
    segment = cell(0); %Initialize Cell
    
    
    Zq = Z;
    Zq(2:end-1,2:end-1) = nan(size(Zq(2:end-1,2:end-1)));
    
    if all(Z(:) > level) || all(Z(:) < level)
        continue; %Skip Level if all Z-Points are above or below level
    elseif ~(all(Zq(~isnan(Zq)) > level) ||  all(Zq(~isnan(Zq)) < level))
        if closeContours
            Z(Zq>level) = level;
        end
    end
    
    
    % Loop through GRID
    for ii = 1:rows-1
        for jj = 1:cols-1
            
            top = ii; bottom = ii+1; left = jj; right = jj+1; %For readability
            %Square vertices
            Zr = Z([top,bottom],[left,right]);
            
            %Skip if all vertices are below or above level
            if all(Zr(:)<level) || all(Zr(:)>level)
                continue;
            end
            
            % Square center and value
            center = [mean([X(top,left) X(bottom,right)]),mean([Y(top,left) Y(bottom,right)])];
            centerMean = mean(Zr(:));
            
            %      vertex 4 +-------------------+ vertex 3
            %               | \               / |
            %               |   \    m=4    /   |
            %               |     \       /     |
            %               |       \   /       |
            %               |  m1   center  m3  |
            %               |       /   \       |
            %               |     /       \     |
            %               |   /    m=2    \   |
            %               | /               \ |
            %      vertex 1 +-------------------+ vertex 2
            
            points = [top left bottom left;... % m1
                bottom left bottom right;...   % m2
                bottom right top right;...     % m3
                top right top left];           % m4
            
            %Loop m1 to m4
            %Each triangle uses 3 points, p1, p2, and the center point
            for m = 1:4 %For all 4 Triangles in rectangle
                section = points(m,:);
                p1 = Z(section(1),section(2)); %Triangle point (1) Value and Coordinates
                p1coords = [X(section(1),section(2)), Y(section(1),section(2))];
                p2 = Z(section(3),section(4)); %Triangle point (2) Value and Coordinates
                p2coords = [X(section(3),section(4)), Y(section(3),section(4))];
                
                %Puts points and center coordinates in a matrix
                coords = [p1coords ; p2coords ; center];
                pointVals = [p1; p2 ; centerMean];
                
                if plotTriangles
                    plot(X,Y,'k.');
                    hold on
                    plot(X([ii,ii+1],[jj,jj+1]),Y([ii,ii+1],[jj,jj+1]),'r.','MarkerSize', 10);
                    plot([coords(:,1) ; coords(1,1)],[coords(:,2) ; coords(1,2)],'r','LineWidth',0.1)
                end
                
                %Points in Triangle on/below/above level
                numOn = sum([p1 p2 centerMean] == level);
                numBelow = sum([p1 p2 centerMean] < level);
                numAbove = sum([p1 p2 centerMean] > level);
                
                if numBelow == 2 &&  numAbove == 1 ||... %case c
                        numBelow == 1 && numAbove == 2   %case f
                    %line is drawn from one edge to another edge of the triangle
                    abVert = coords([p1 p2 centerMean] > level,:);
                    abVal = pointVals([p1 p2 centerMean] > level,:);
                    belVert = coords([p1 p2 centerMean] < level,:);
                    belVal = pointVals([p1 p2 centerMean] < level,:);
                    
                    if numAbove == 1
                        %Interpolate between point above and points below
                        X1 = (level-abVal(1))*(belVert(1,1)-abVert(1,1))/(belVal(1)-abVal(1))+abVert(1,1);
                        Y1 = (level-abVal(1))*(belVert(1,2)-abVert(1,2))/(belVal(1)-abVal(1))+abVert(1,2);
                        %Interpolate
                        X2 = (level-abVal(1))*(belVert(2,1)-abVert(1,1))/(belVal(2)-abVal(1))+abVert(1,1);
                        Y2 = (level-abVal(1))*(belVert(2,2)-abVert(1,2))/(belVal(2)-abVal(1))+abVert(1,2);
                        
                    else %elseif numBelow == 1
                        %Interpolate between point below and points above
                        X1 = (level-belVal(1))*(abVert(1,1)-belVert(1,1))/(abVal(1)-belVal(1))+belVert(1,1);
                        Y1 = (level-belVal(1))*(abVert(1,2)-belVert(1,2))/(abVal(1)-belVal(1))+belVert(1,2);
                        %Interpolate
                        X2 = (level-belVal(1))*(abVert(2,1)-belVert(1,1))/(abVal(2)-belVal(1))+belVert(1,1);
                        Y2 = (level-belVal(1))*(abVert(2,2)-belVert(1,2))/(abVal(2)-belVal(1))+belVert(1,2);
                        
                    end
                    segment{end+1} = round([X1 Y1 ; X2 Y2],4); %Add line to cell array
                    
                elseif (numBelow == 1 && numOn == 2) ||... %case d
                        (numOn == 2 && numAbove == 1)      %case h
                    %line is drawn between the two vertices that lie on the contour level
                    segment{end+1} = round(coords([p1 p2 centerMean] == level,:),4); %Add line to cell array
                elseif numBelow == 1 && numAbove ==1 && numOn == 1 %case e
                    %line be drawn from the vertex on the contour level to a point on the opposite edge
                    onVert = coords([p1 p2 centerMean] == level,:);
                    oppVert = coords([p1 p2 centerMean] ~= level,:);
                    oppVal =  pointVals([p1 p2 centerMean] ~= level);
                    %Interpolate between point above and point below level
                    vertOppX = (level-oppVal(1))*(oppVert(2,1)-oppVert(1,1))/(oppVal(2)-oppVal(1))+oppVert(1,1);
                    vertOppY = (level-oppVal(1))*(oppVert(2,2)-oppVert(1,2))/(oppVal(2)-oppVal(1))+oppVert(1,2);
                    %Add vertex on level and interpolated point
                    segment{end+1} = round([onVert;[vertOppX, vertOppY]],4); %Add line to cell array
                else
                    continue;
                end
                %Deletes 'lines' that have 0 length due to rounding
                if all(segment{end}(1,:) == segment{end}(2,:))
                    segment(end) = [];
                end
            end
        end
    end
    
    if isempty(segment)
        continue;
    end
    
    
    %% SECTION 2: DELETE LINE SEGMENT COPIES
    %trackr is an array that tracks the status of the "segment" cell array
    %throughout the sorting progress.  Line segments are represented by cell
    %elements in "segment" cell array.  In this section, we change the
    %duplicate segment to NaNs so they cannot be matched with anything later,
    %meanwhile recording that segments index in trackr.  They are then deleted.
    trackr = zeros(length(segment),1); %initialize
    for ii = 1:length(segment)
        for jj = 1:length(segment)
            if ii == jj
                continue;
            end
            %Flip Line Segment
            opp = flip(segment{jj});
            %See if the segment matches another exactly or its reverse
            if all(segment{ii}(:) == segment{jj}(:)) || all(segment{ii}(:) == opp(:))
                segment{jj} = [NaN NaN; NaN NaN];
                trackr(jj) = 1;
            end
        end
    end
    segment = segment(~logical(trackr)); %Remove used lines
    trackr = trackr(~logical(trackr)); %Shorten line tracker so it still works later
    
    %% SECTION 3: SORT LINES INTO DISTINCT CONTOURS
    
    cN = 0;
    while ~isempty(segment) %Every line must belong to a contour DUH
        SkipIndex = 1;
        cN = cN+1; % Cell number +1
        contourCell{cN} = [];
        coordinate = segment{1}(1,:);
        match = false;
        while ~match
            for ii = 1:length(segment)%Iterate through all lines
                if SkipIndex == ii || isempty(segment{ii})
                    %Do Nothing
                elseif all(coordinate == segment{ii}(1,:)) % Check first vertex in line
                    coordinate = segment{ii}(2,:); % Make other coordinate in the line the next one to compare
                    contourCell{cN}(end+1,:) = segment{ii}(1,:); % Add matched coordinate to contour
                    segment{ii} = [];
                    trackr(ii) = cN; % Modify list
                    SkipIndex = ii; % Make this index the one to skip in the next iteration
                    break;
                elseif all(coordinate == segment{ii}(2,:)) % Check second vertex in line
                    coordinate = segment{ii}(1,:);
                    contourCell{cN}(end+1,:) = segment{ii}(2,:);
                    segment{ii} = [];
                    trackr(ii) = cN;
                    SkipIndex = ii;
                    break;
                end % if all(vert
                
                if ii == length(segment) %No matches signaling end of the contour
                    %Record last coordinate
                    contourCell{cN}(end+1,:) = coordinate;
                    %Throw copy first point at the end and declare no match
                    if isempty(segment{1})
                        %Add to contour array
                        C = [C, [level ; length(contourCell{cN})], contourCell{cN}'];
                        match = true;
                    else
                        contourCell{cN} = flip(contourCell{cN},1); %flip on 1st dim
                        coordinate = segment{1}(2,:);
                        segment{1} = [];
                        trackr(1) = cN;
                        SkipIndex = 1;
                    end
                end
            end % for 1:length(segment)
            
            if plotFlag
                hold on;
                for qq = 1:length(contourCell)
                    plot(contourCell{qq}(:,1),contourCell{qq}(:,2),'r','LineWidth',2)
                end
                drawnow
            end
            
        end % while ~match
        segment = segment(~logical(trackr)); %Remove used lines
        trackr = trackr(~logical(trackr)); %Shorten line tracker so it still works later
    end %while ~all(trackr)
    
    if ~plotFlag
        hold on;
        for qq = 1:length(contourCell)
            plot(contourCell{qq}(:,1),contourCell{qq}(:,2),'r','LineWidth',1)
        end
        drawnow
    end
    
end
end
