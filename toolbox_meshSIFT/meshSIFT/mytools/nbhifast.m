function neighborhood = nbhifast(data,ncols,index,d,zmin,dpr)
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
% dpr = 0 means the datamatrix is arranged per column, 1 means per row.
% neighborhood = a list of 3D-points in the neighborhood of point p.
% Background points will not be included.

%% Initialisation
if nargin < 6
    dpr = 0;
    if nargin < 5
        zmin = -1E9;
        if nargin < 4
            d = 3;
            if nargin < 3
                disp('you did not give the required number of parameters.')
            end
        end
    end
end
[len,nc] = size(data);
if nc>3
    data = data(:,1:3);
end
nrows = len/ncols;
if dpr == 0
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
nel = (d1+d4+1)*(d2+d3+1);
temp = zmin*ones(nel,3); 
i = index - d2 - d1*ncols;
k = 1;
while i <= (index + d3 + d4*ncols)
    j = i;
    while j <= i+d2+d3
        if data(j,3) ~= zmin
            temp(k,:) = data(j,:);
            k = k+1;
        end
        j = j + 1;
    end
    i = i + ncols;
end
neighborhood = temp(1:k-1,:);
end