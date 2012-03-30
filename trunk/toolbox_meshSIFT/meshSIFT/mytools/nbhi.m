function neighborhood = nbhi(data,ncols,index,d,zmin,dpc)
% generates the neighborhood of a certain point (index given) 
% in a depth scan, eliminating background points
% Author: Chris Maes
% 2009/09

% Usage
% data = matrix with face data (array of 3D data points)
% index = index of the point of intrest (row nr)
% d = distance to go from point p. d=3 will generate a 7 by 7
%     neighbor region (d = 3 by default)
% zmin = z-value indicating background (-1E9 by default)
% dpc = 0 means the datamatrix is arranged per row. 1 means per column.
% neighborhood = a list of 3D-points in the neighborhood of point p.
% Background points will not be included.

%% Initialisation
if nargin < 6
    dpc = 0;
    if nargin < 5
        zmin = -1E9;
        if nargin < 4
            d = 3;
            if nargin < 3
                warning('MATLAB:ParamAmbiguous','you did not give the required number of parameters.')
            end
        end
    end
end
len = size(data,1);
nrows = len/ncols;
if dpc == 1
    temp = nrows;
    nrows = ncols;
    ncols = temp;
end
if round(nrows) ~= nrows
    disp('The dimensions of the datamatrix and number of columns are not correct.')
end
if index<=0 || index > len
    disp('the index is not within the matrix.')
end

%% Calculate distances to go from index

rowindex = floor((index-1)/ncols)+1; % index of the row of the point of intrest
restc = rem(index,ncols); % position in that row
if restc == 0
    restc = ncols;
end
if rowindex > d
    d1 = d; % number of positions to go up in the matrix
else
    d1 = rowindex - 1;
end
if restc > d
    d2 = d; % number of positions to go to the left in the matrix
else
    d2 = restc - 1;
end
if restc + d <= ncols
    d3 = d; % number of positions to go to the left in the matrix
else
    d3 = ncols - restc;
end
if rowindex + d <= nrows
    d4 = d; % number of positions to go down in the matrix
else
    d4 = nrows - rowindex;
end

%% Put the points of intrest in a matrix
temp = []; 
i = index - d2 - d1*ncols;
while i <= (index + d3 + d4*ncols)
    j = i;
    while j <= i+d2+d3
        if data(j,1) ~= zmin
            temp = [temp; data(j,:)]; %#ok<AGROW>
        end
        j = j + 1;
    end
    i = i + ncols;
end
neighborhood = temp;
end