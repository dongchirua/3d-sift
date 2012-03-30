function indices = nbhmi(vertex,faces,index,r)

% generates the neigborhood of a point at a given index in the vertex
% array.
% vertex = nx3 matrix with the vertex points from a mesh
% faces = mx3 matrix with the corner indices of the faces of a mesh
% index = index of point of intrest in vertex
% r = the radius of the neighborhood
% nbh = list of indices in the neighborhood of vertex(index)
% Author: Chris Maes
% 2009/10

%% calculation indices
D = perform_fast_marching_mesh(vertex, faces, index);
indices = find(D<=r);
end

    