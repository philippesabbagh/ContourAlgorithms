%%%%%%%%% CONTOUR TEST MATLAB Edition %%%%%%%%%
%%%%%%%%% CONTOUR TEST 1 %%%%%%%%%
close all; clear;
[X,Y] = meshgrid(-1.5:0.1:1.5 , -1.5:0.1:1.5);
Z = 1./((X.^2+(Y-0.842).*(Y + 0.842)).^2+(X.*(Y-0.842)+X.*(Y-0.842)).^2);
LEVELS = [0.2 1 1.5 2 3];

C = contourPhilippeMATLAB(X,Y,Z,LEVELS);

if length(LEVELS) > 1
    contour(X,Y,Z,LEVELS, 'b.-', 'LineWidth',1);
else
    contour(X,Y,Z,[LEVELS LEVELS], 'b.-', 'LineWidth',1);
end

%%%%%%%%% CONTOUR TEST 2 %%%%%%%%%
close all; clear;
[X,Y] = meshgrid(-1.5:0.1:1.5 , -1.5:0.1:1.5);
Z = sin((X + Y.^2));
LEVELS = [-1:0.2:1];

C = contourPhilippeMATLAB(X,Y,Z,LEVELS);

if length(LEVELS) > 1
    contour(X,Y,Z,LEVELS, 'b-.', 'LineWidth',1);
else
    contour(X,Y,Z,[LEVELS LEVELS], 'b-.', 'LineWidth',1);
end

%%%%%%%%% CONTOUR TEST 3 %%%%%%%%%
close all; clear;
[X,Y] = meshgrid(-1.5:0.5:1.5 , -1.5:0.5:1.5);

Z = [...
    0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 1, 1;...
    1, 1, 1, 1, 1, 1, 0;...
    0, 0, 1, 0, 0, 1, 0;...
    0, 0, 1, 0, 0, 1, 0;...
    0, 0, 1, 1, 1, 1, 0;...
    0, 0, 1, 0, 0, 0, 0 ...
    ];

LEVELS = [0.2 0.5 0.7 0.8];

C = contourPhilippeMATLAB(X,flip(Y),Z,LEVELS);

if length(LEVELS) > 1
    contour(X,flip(Y),Z,LEVELS, 'b.-', 'LineWidth',1);
else
    contour(X,flip(Y),Z,[LEVELS LEVELS], 'b.-', 'LineWidth',1);
end

%%%%%%%%% CONTOUR TEST 4 %%%%%%%%%
close all; clear;
[X,Y] = meshgrid(-1.5:0.5:1.5 , -1.5:0.5:1.5);

        Z = [...
            0, 1, 0, 1, 1, 0, 0;...
            1, 1, 0, 1, 0, 0, 1;...
            0, 1, 0, 0, 0, 1, 0;...
            0, 1, 0, 1, 1, 0, 0;...
            0, 1, 0, 1, 1, 0, 0;...
            0, 1, 0, 0, 0, 0, 0;...
            0, 0, 0, 0, 0, 0, 1 ...
            ];

LEVELS = [0.25 0.6 0.9];

C = contourPhilippeMATLAB(X,flip(Y),Z,LEVELS);

if length(LEVELS) > 1
    contour(X,flip(Y),Z,LEVELS, 'b.-', 'LineWidth',1);
else
    contour(X,flip(Y),Z,[LEVELS LEVELS], 'b.-', 'LineWidth',1);
end

%%%%%%%%% CONTOUR TEST 5 %%%%%%%%%
close all; clear;
[X,Y] = meshgrid(0:1 , 0:1);

        Z = [...
            1, 0.4; ...
            0, 0.8 ...
            ];

LEVELS = [0.6];

C = contourPhilippeMATLAB(X,flip(Y),Z,LEVELS);

if length(LEVELS) > 1
    contour(X,flip(Y),Z,LEVELS, 'b.-', 'LineWidth',1);
else
    contour(X,flip(Y),Z,[LEVELS LEVELS], 'b.-', 'LineWidth',1);
end
