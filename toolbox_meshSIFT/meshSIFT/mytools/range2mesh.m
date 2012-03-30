function [vertex,faces] = range2mesh(data,nrows,ncols,zmin,nvert,df)
% [vertex,faces] = range2mesh(data,nrows,ncols,zmin,nvert,df)
%
% This function converts a bosphorus range-image into a mesh, ignoring all
% background points with value zmin.
% nvert = number of vertices you would like to retain after downsampling, 
% df = downsamplefactor;
% ignoring zmin values
% 
% Author: Chris Maes, Dirk Smeets
% 2009/10

%% Initialisation
if nargin < 5
    nvert = size(data,1);
    if nargin < 4
        warning('MATLAB:paramAmbiguous','you did not give the required number of parameters.')
    end
end
if nvert > size(data,1)
    nvert = size(data,1);
end
dc = size(data,2);
if dc>3
    data = data(:,1:3);
elseif dc<3
    warning('MATLAB:paramAmbiguous','The data matrix must contain 3 columns')
end   
npoints = nrows*ncols;
if nargin >5
    iv = df;
    nvert = 1;
else
    iv = floor(sqrt(npoints/nvert)); % interval
end
nvertices = 0;
while nvertices < nvert && iv > 0
    startc = ceil(rem(ncols,iv)/2);
    endc = floor(rem(ncols,iv)/2);
    if startc == 0
        startc = ceil(iv/2);
        endc = iv-startc;
    end
    startr = ceil(rem(nrows,iv)/2);
    endr = floor(rem(nrows,iv)/2);
    if startr == 0
        startr = ceil(iv/2);
        endr = iv-startr;
    end
    mpr = ceil(2*ncols/iv+6); % maximum number of faces per row (+6)

%% Calculation
    vertex = zeros(npoints,3);
    j = 1;
    faces = zeros(3*npoints,3);
    k = 1;
    for r = startr:iv:nrows-iv
        for c = startc:iv:ncols
            if c == ncols-endc
                index = (r-1)*ncols+c;
                if data(index,3) ~= zmin
                    vertex(j,:) = data(index,:);
                    j = j + 1;
                    if k < mpr+2
                        faces(faces==index)=j-1;
                    else
                        ind = false(k-1,3);
                        ind(k-mpr-1:k-1,:) = faces(k-mpr-1:k-1,:)==index;
                        faces(ind(:,1)) = j-1;
                        faces(ind(:,2),2) = j-1;
                        faces(ind(:,3),3) = j-1;
                    end
                end
            else
                index = (r-1)*ncols+c;
                if data(index,3) ~= zmin
                    vertex(j,:) = data(index,:);
                    j = j + 1;
                    if data(index+iv*ncols+iv,3) ~= zmin
                        if data(index+iv,3) ~= zmin
                            faces(k,:) = [index index+iv*ncols+iv index+iv];
                            k = k + 1;
                        end
                        if data(index+iv*ncols,3) ~= zmin
                            faces(k,:) = [index index+iv*ncols index+iv*ncols+iv];
                            k = k + 1;
                        end
                    else
                        if data(index+iv,3) ~= zmin
                            if data(index+iv*ncols,3) ~= zmin
                                faces(k,:) = [index index+iv*ncols index+iv];
                                k = k + 1;
                            end
                        end
                    end
                    if k < mpr+2
                        faces(faces==index)=j-1;
                    else
                        ind = false(k-1,3);
                        ind(k-mpr-1:k-1,:) = faces(k-mpr-1:k-1,:)==index;
                        faces(ind(:,1)) = j-1;
                        faces(ind(:,2),2) = j-1;
                        faces(ind(:,3),3) = j-1;
                    end
                else
                    if data(index+iv*ncols,3) ~= zmin
                        if data(index+iv,3) ~= zmin
                            if data(index+iv*ncols+iv,3) ~= zmin
                                faces(k,:) = [index+iv index+iv*ncols index+iv*ncols+iv];
                                k = k + 1;
                            end
                        end
                    end
                end
            end
        end
    end
    r = nrows-endr;
    for c = startc:iv:ncols
        index = (r-1)*ncols+c;
        if data(index,3) ~= zmin
            vertex(j,:) = data(index,:);
            j = j + 1;
            ind = false(k-1,3);
            if k < mpr+2
                faces(faces==index)=j-1;
            else
                ind(k-mpr-1:k-1,:) = faces(k-mpr-1:k-1,:)==index;
                faces(ind(:,1)) = j-1;
                faces(ind(:,2),2) = j-1;
                faces(ind(:,3),3) = j-1;                
            end
        end
    end
    vertex = vertex(1:j-1,:);
    faces = faces(1:k-1,:);
    nvertices = j-1;
    iv = iv-1;
end
end