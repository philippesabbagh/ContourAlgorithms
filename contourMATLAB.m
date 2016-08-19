% Philippe Sabbagh
% Draper

function [C] = contourPhilippeMATLAB(X,Y,Z,LEVELS)
plotFlag = 0;
closeContours = 0; %EXPERIMENTAL
C = [];
[rows, cols] = size(X);
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
    
    
    % Loop through GRID
    for ii = 1:rows-1
        for jj = 1:cols-1
            
            tp = ii; bm = ii+1; lt = jj; rt = jj+1; %For readability
            %Square vertices
            Zr = Z([tp,bm],[lt,rt]);
            %Skip if all vertices are below or above level
            if all(Zr(:)<level) || all(Zr(:)>=level)
                continue;
            end
            
            %               +-p2-------( M4 )-------p1-+
            %               |                          |
            %               p1                         p2
            %               |                          |
            %               |                          |
            %             ( M1 )                     ( M3 )
            %               |                          |
            %               |                          |
            %               p2                         p1
            %               |                          |
            %               +-p1-------( M2 )-------p2-+ 
            
            points = [tp lt bm lt;... % m1
                      bm lt bm rt;... % m2
                      bm rt tp rt;... % m3
                      tp rt tp lt];   % m4
                  
%             plot(X([tp,bm],[lt,rt]),Y([tp,bm],[lt,rt]),'r.')

            vert = [];
            for m = 1:4 %Loop m1 to m4
                section = points(m,:);
                p1 = Z(section(1),section(2));
                p1coords = [X(section(1),section(2)), Y(section(1),section(2))];
                p2 = Z(section(3),section(4));
                p2coords = [X(section(3),section(4)), Y(section(3),section(4))];
                
                if abs(p2-p1) > 10e4
                    disp(['Gradient High: Rounding error may occur. Row: ' ...
                        num2str(ii) ' Col: ' num2str(jj)])
                end
                
                if (p1>level && p2>level) ||...
                        (p1<level && p2<level) ||...
                        (p1==level && p2==level)
                    continue;
                elseif p1 == Inf || p2 == Inf
                    if p1 == Inf
                        vert = [vert; round(p2coords,4)];                                             
                    else
                        vert = [vert; round(p1coords,4)];   
                    end
                    
                elseif p1 == -Inf || p2 == -Inf
                    if p1 == -Inf                        
                        vert = [vert; round(p2coords,4)];                                             
                    else
                        vert = [vert; round(p1coords,4)];   
                    end
                else
                    vertX = round((level-p1)*(p2coords(1)-p1coords(1))/(p2-p1)+p1coords(1),4);
                    vertY = round((level-p1)*(p2coords(2)-p1coords(2))/(p2-p1)+p1coords(2),4);
                    vert = [vert; vertX vertY];
                end
                
            end
            
            
            %Make decisions based on number of vertices
            [numVert,~] = size(vert);
            numUnique = length(unique(vert,'rows'));
            
            if numVert == 1
                continue;
            elseif numVert == 2
                segment{end+1} = vert; %Add line to cell array
            elseif numVert == 3
                if numUnique == 2
                    vert = unique(vert,'rows'); %Delete duplicate vertices
                    segment{end+1} = vert; %Add line to cell array
                elseif numUnique == 3
                    error('Investigate 3 vertex case')
                end
                
            elseif numVert == 4
                if numUnique == 1
                    continue; %Shouldn't Happen Really
                elseif numUnique == 2
                    vert = unique(vert,'rows'); %Delete duplicate vertices
                    segment{end+1} = vert; %Add line to cell array
                    continue
                elseif numUnique == 3
                    error('Investigate 3 vertex case')
                else
                    distance = [0;0;0;0];
                    for vv = 2:4
                        distance(vv) = ((vert(vv,1)-vert(1,1))^2+(vert(vv,2)-vert(1,2))^2)^0.5;
                    end
                    diagDistance = distance([2,4]);
                    %Connect to closest vertex
                    if diagDistance(1) < diagDistance(2)
                        segment{end+1} = [vert(1,:) ; vert(2,:)]; %Add segment
                        vert([1,2],:) = []; %Delete used points
                    elseif diagDistance(2) < diagDistance(1)
                        segment{end+1} = [vert(1,:) ; vert(4,:)]; %Add segment
                        vert([1,4],:) = []; %Delete used points
                    else %If side vertices are equally far away.
                        if Z(ii,jj) > level
                            segment{end+1} = [vert(1,:) ; vert(2,:)]; %Add segment
                            vert([1,2],:) = []; %Delete used points
                        else
                            segment{end+1} = [vert(1,:) ; vert(4,:)]; %Add segment
                            vert([1,4],:) = []; %Delete used points
                        end
                    end
                    segment{end+1} = vert; %Save Remaining points as segment;
                end
                
            else
                error('Investigate weird case :(')
            end
            
        end
    end
    
    if isempty(segment)
        continue;
    end
    
    
    %% SECTION 3: SORT LINES INTO DISTINCT CONTOURS
    trackr = zeros(length(segment),1); %initialize
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
