% Philippe Sabbagh
% Draper

function [C] = contourPhilippeMATLAB(X,Y,Z,LEVELS)
plotFlag = 0;
closeContours = 0; %EXPERIMENTAL
C = [];
[rows, cols] = size(X);
% plot(X,Y,'k.');
% Y = flip(Y);
figure; grid on; axis([min(X(:)) max(X(:)) min(Y(:)) max(Y(:))]); hold on;
%% SECTION 1: PHIL SUBROUTINE
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
    
    vert = [];
    
    % Loop through GRID
    for ii = 1:rows-1
        for jj = 1:cols-1
            
            top = ii; bottom = ii+1; left = jj; right = jj+1; %For readability
            %Square vertices
            Zr = Z([top,bottom],[left,right]);
            
            %Skip if all vertices are below or above level
            if all(Zr(:)<=level) || all(Zr(:)>=level)
                continue;
            end            
            
            %      vertex 4 +--------m4----------+ vertex 3
            %               |                   |
            %               |                   |
            %               |                   |
            %               |                   |
            %              m1                   m3
            %               |                   |
            %               |                   |
            %               |                   |
            %               |                   |
            %      vertex 1 +--------m2---------+ vertex 2
            
            points = [top left bottom left;... % m1
                bottom left bottom right;...   % m2
                bottom right top right;...     % m3
                top right top left];           % m4
            
            %Loop m1 to m4
            %Each triangle uses 3 points, p1, p2, and the center point
            
            vert = [];
            
            for m = 1:4 %For all 4 Triangles in rectangle
                section = points(m,:);
                p1 = Z(section(1),section(2));
                p1coords = [X(section(1),section(2)), Y(section(1),section(2))];
                p2 = Z(section(3),section(4));
                p2coords = [X(section(3),section(4)), Y(section(3),section(4))];
                
                %Puts points and center coordinates in a matrix
                coords = [p1coords ; p2coords];
                pointVals = [p1; p2];
                
                if (p1>level && p2>level) ||...
                        (p1<level && p2<level) ||...
                        (p1==level && p2==level)
                    continue;
                else
                    vertX = round((level-p1)*(p2coords(1)-p1coords(1))/(p2-p1)+p1coords(1),4);
                    vertY = round((level-p1)*(p2coords(2)-p1coords(2))/(p2-p1)+p1coords(2),4);
                    vert = [vert; vertX vertY];
                end    
                
            end
            
            
            %Make decisions based on number of vertices
            vert = unique(vert,'rows'); %Delete duplicate vertices
            [numVert,~] = size(vert);
            
            if numVert == 1
                continue;
            elseif numVert == 2 %Happy Happy
                segment{end+1} = vert; %Add line to cell array                
            elseif numVert == 3
                error('Investigate 3 vertex case')
            elseif numVert == 4
                distance = [0;0;0;0];
                for vv = 2:4
                    distance(vv) = ((vert(vv,1)-vert(1,1))^2+(vert(vv,2)-vert(1,2))^2)^0.5;                    
                end
                [~,II] = min(distance(2:4));
                segment{end+1} = [vert(1,:) ; vert(II+1,:)]; %Add segment 
                
                
                vert([1,II+1],:) = []; %Delete used points
                segment{end+1} = vert; %Save Remaining points as segment;
                
            else
                error('Investigate weird case :(')
            end
            
        end
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
