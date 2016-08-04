% Desperately trying to make contours
% Philippe Sabbagh
%
% CONREC subroutine by Paul Bourke July 1987 http://paulbourke.net/papers/conrec/

function [C] = contourPrototype(X,Y,Z,LEVELS)
[cols, rows] = size(X);
% plot(X,Y,'k.');
% Y = flip(Y);
figure; grid on; axis equal; hold on;
%% CONREC SUBROUTINE
for level = LEVELS
    vert = cell(0); %Initialize Cell
    % Loop through GRID
    for ii = 1:rows-1
        for jj = 1:cols-1
            top = ii; bottom = ii+1; left = jj; right = jj+1; %For readability
            %Square vertices
            Zr = Z([ii,ii+1],[jj,jj+1]);
            % Square center and value
            center = [mean([X(ii,jj) X(ii+1,jj+1)]) , mean([Y(ii,jj) Y(ii+1,jj+1)])];
            centerMean = mean(Zr(:));
            %Skip if all vertices are below or above level
            if all(Zr(:)<level) || all(Zr(:)>level)
                continue;
            end
            
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
                    vert{end+1} = round([X1 Y1 ; X2 Y2],4); %Add line to cell array
                    
                elseif (numBelow == 1 && numOn == 2) ||... %case d
                        (numOn == 2 && numAbove == 1)      %case h
                    %line is drawn between the two vertices that lie on the contour level
                    vert{end+1} = round(coords([p1 p2 centerMean] == level,:),4); %Add line to cell array
                elseif numBelow == 1 && numAbove ==1 && numOn == 1 %case e
                    %line be drawn from the vertex on the contour level to a point on the opposite edge
                    onVert = coords([p1 p2 centerMean] == level,:);
                    oppVert = coords([p1 p2 centerMean] ~= level,:);
                    oppVal =  pointVals([p1 p2 centerMean] ~= level);
                    %Interpolate between point above and point below level
                    vertOppX = (level-oppVal(1))*(oppVert(2,1)-oppVert(1,1))/(oppVal(2)-oppVal(1))+oppVert(1,1);
                    vertOppY = (level-oppVal(1))*(oppVert(2,2)-oppVert(1,2))/(oppVal(2)-oppVal(1))+oppVert(1,2);
                    %Add vertex on level and interpolated point
                    vert{end+1} = round([onVert;[vertOppX, vertOppY]],4); %Add line to cell array
                else
                    % Gaps in their contour lines, this should of
                    % course never happen. There is however a pathological
                    % case that all local contouring algorithms suffer from
                    % (local meaning that they only use information in the immediate
                    % vicinity to determine the contour lines). The problem arises when
                    % all four vertices of a grid cell have the same value as the contour
                    % level under consideration. There are a number of strategies that
                    % can be employed to overcome this special event, the correct way
                    % is to consider a larger region in order to join up the contours
                    % on either side of the problem cell. CONREC doesn't do this
                    % and just leaves the cell without any contour lines thus resulting
                    % in a gap. This special case essentially never happens for real
                    % values data, it is most commonly associated with integer height
                    % datasets. The simplest solution is to offset the contour levels
                    % being drawn by a very small amount.
                    continue;
                end                
            end            
        end
    end
    
    if isempty(vert)
        warning(['Level ' num2str(level) ' does not produce a contour'])
        continue;
    end  
    
    %% SORT LINES INTO DISTINCT CONTOURS
    C = [];
    cN = 0;
    ind = 1;
    vl = zeros(length(vert),1); %initialize
    while ~all(vl) %Every line must belong to a contour DUH
        cN = cN+1; % Cell number +1
        contourCell{cN} = []; %New Cell for a new contour
        %Choose first vertex in first line as arbitrary start point
        coordinate = vert{1}(1,:);
        FC = vert{ind}(1,:);
        match = true;
        while match
            for ii = 1:length(vert)%Iterate through all lines                
                if ind == ii
                    continue;
                elseif all(coordinate == vert{ii}(1,:)) % Check first vertex in line
                    coordinate = vert{ii}(2,:); % Make other coordinate in the line the next one to compare
                    contourCell{cN}(end+1,:) = vert{ii}(1,:); % Add matched coordinate to contour
                    vert{ii}(1,:) = [];  % Make into bread so it isn't used again
                    vl(ii) = cN; % Modify list
                    ind = ii; % Make this index the one to skip in the next iteration
                    break;
                elseif all(coordinate == vert{ii}(2,:)) % Check second vertex in line
                    coordinate = vert{ii}(1,:);
                    contourCell{cN}(end+1,:) = vert{ii}(2,:);
                    vert{ii}(2,:) = [NaN NaN];
                    vl(ii) = cN;
                    ind = ii;
                    break;
                end % if all(vert
                
                if ii == length(vert) %No matches signaling end of the contour
                    %Throw copy first point at the end and declare no match
                    contourCell{cN}(end+1,:) = vert{1}(1,:);
                    match = false;
                    %Add to contour array
                    C = [C, [level ; length(contourCell{cN})], contourCell{cN}'];
                end
            end % for 1:length(vert)
            
        end % while match
        %delete empty cell arrays
        logicarray = ~cellfun(@isempty, vert);
        vert = vert(logicarray);
%         vert = vert(~logical(vl)); %Remove used lines
%         vl = vl(~logical(vl)); %Shorten line tracker so it still works later
    end %while ~all(vl)
    
    for ii = 1:length(contourCell)
        plot(contourCell{ii}(:,1),contourCell{ii}(:,2),'k')
    end
end
end

% Contour Test
close all; clear;
[X,Y] = meshgrid(-1.5:0.1:1.5 , -1.5:0.1:1.5);
Z = 1./((X.^2+(Y-0.842).*(Y + 0.842)).^2+(X.*(Y-0.842)+X.*(Y-0.842)).^2);
LEVELS = 1;

C = contourPrototype(X,Y,Z,LEVELS);
% if length(LEVELS) > 1
%     contour(X,Y,Z,LEVELS, 'r');
% else
%     contour(X,Y,Z,[LEVELS LEVELS], 'r');
% end

% legend('Phil Contour','Matlab built in contour.m')
