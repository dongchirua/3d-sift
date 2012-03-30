function nbh = nbhm(vertex,faces,index,r)

% generates the neigborhood of a point at a given index in the vertex
% array.
% vertex = nx3 matrix with the vertex points from a mesh
% faces = mx3 matrix with the corner indices of the faces of a mesh
% index = index of point of intrest in vertex
% r = the radius of the neighborhood (geodetic distance)
% nbh = list of points in the neighborhood of vertex(index)
% Author: Chris Maes
% 2009/10
%% initialisation
if nargin<4
    r = 1;
    if nargin < 3
        warning('MATLAB:ParamAmbiguous','you did not give the required number of parameters.')
    end
end

%% calculation indices
D = perform_fast_marching_mesh(vertex, faces, index);
indices = find(D<=r);

%% convert to 3D points
len = length(indices);
nbh = zeros(len,3);
j = 1;
for i = 1:len
    nbh(j,:)= vertex(indices(i),:);
    j = j+1;
end
end

    