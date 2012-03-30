function neighborhood = nbh(A,pr,pc,d,zmin)
% generates the neighborhood of a certain point in a depth scan
% eliminating background points
% Author: Chris Maes
% 2009/09

% Usage
% A = matrix with face data
% pr and pc are the row- and column-indices of the point of intrest
% d = distance to go from point p. d=3 will generate a 7 by 7
% neighbor region (d = 3 by default)
% zmin = z-value indicating background (-1E9 by default)
% neighborhood = a list of 3D-points in the neighborhood of point p.
% Background points will not be included.

%% Initialisation
if nargin == 4
    zmin = -1E9;
elseif nargin == 3
    d = 3;
    zmin = -1E9;
else
    disp('you did not give the required number of parameters')
end
[nrows,ncols] = size(A);
startr = pr - d; % start position of neighborhood in x-direction
startc = pc - d; % start position of neighborhood in y-direction
stopr = pr + d; % stop position of neighborhood in x-direction
stopc = pc + d; % stop position of neighborhood in y-direction

%% fix start and stop for points near the edge
if startr < 1
    startr = 1;
end
if stopr > nrows
    stopr = nrows;
end
if startc < 1
    startc = 1;
end
if stopc > ncols
    stopc = ncols;
end

%% Put the points of intrest in a matrix
neighborhood = []; 
for i = startr:stopr
    for j = startc:stopc
        if A(i,j) ~= zmin
            neighborhood = [neighborhood; i j A(i,j)]; %#ok<AGROW>
        end
    end
end
end