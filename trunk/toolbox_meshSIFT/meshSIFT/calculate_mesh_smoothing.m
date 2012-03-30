function vertices = calculate_mesh_smoothing(vertex,faces,nscales,prior_smoothing)
% calculate mesh smoothing for different scales
% 
% INPUTS
% ------
% Vertex = Mx3 matrix,
% face = Nx3 matrix
% nscales = Qx1 vector containing number of smoothings steps for the
%           desired result
% prior_smoothing = number of smoothing steps before calculating the
%                   smoothed versions of the mesh
%                 = 0 by default
%
% OUTPUTS
% -------
% vertices = Mx3xQ matrix containing vertex for every mesh on the desired
%            scale in [nscales]
%
% author: Chris Maes
% 2010/04

%% initialisation
if nargin < 4
    prior_smoothing = 0;
    if nargin < 3
        error('3 input arguments needed')
    end
end
nvert = length(vertex);
W = triangulation2adjacency(faces) + speye(nvert); % will be used as a smoothing filter
D = spdiags(full(sum(W,2).^(-1)),0,nvert,nvert);
W = D*W; % normalise rows of filter

temp = vertex;
%% prior smoothing
for i = 1:prior_smoothing
    temp = W*temp;
end

%% calculate smoothed meshes
vertices = zeros(nvert,3,numel(nscales));
nscale = 0;
for i = 1:numel(nscales)
    while nscale < nscales(i)
        nscale = nscale+1;
        temp = W*temp; % smooth with smoothing filter [nscales] times
    end
    vertices(:,:,i) = temp; % save result
end
end