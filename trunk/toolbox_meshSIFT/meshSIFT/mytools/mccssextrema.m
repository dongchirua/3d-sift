function mccssextrema(meshpath,savepath,type)
% give a mesh (vertex and faces), this function returns scale-space extrema
if nargin < 3
    type = 'H'
    if nargin < 2
        error 'apropriate inputs needed'
    end
end
load(meshpath) % this must contain vertex and faces as described below

% INPUTS
% ------
% Vertex = Mx3 matrix,
% faces = Nx3 matrix

% params (optional) = structure containing optional parameters:
% params.k = number of scales per doubling
%          = 4 by default
% params.type = type of extrema to use (K,H,C,S,KK,D)
%             = 'H' by default  
% params.curvature_smoothing = smoothing used for calculation of curvatures
%                            = 3 by default
% params.startscale = number of smoothing steps to start scale space from
%                   = 5 by default
% params.nscalesteps = number of scales where extrema are to be detected,
%                       starting from [startscale] smoothing steps.
%                    = 5 by default (so 8 smoothing scales will be calculated)

% OUTPUT
% ------
% extr = Qx2 matrix. Every row contains the index of the scale-space
% extremum and the assigned scale: [index scale]
% scales = 1xS array containing the scales of smoothing for the mesh (so the
% scales where scale-space extrema were sought + 3 extra scales)
% vertices = Mx3xS matrix, containing vertex for every scale in [scales].

% Author: Chris Maes
% 2010/04

nvert = size(vertex,1);
%% parameters
params = [];
k = getoptions(params,'k',4);
% type = upper(getoptions(params,'type','H'));
options.curvature_smoothing = getoptions(params,'curvature_smoothing',3);
nscalesteps = getoptions(params,'nscalesteps',5);
startscale = getoptions(params,'startscale',5);
prior_smoothing = getoptions(params,'prior_smoothing',1);

%% calculation of scales
W = triangulation2adjacency(faces);
idx = find(triu(W,1));
rows = rem(idx,size(W,1));
cols = ceil(idx/size(W,1));
dsts = sqrt(sum((vertex(rows,:)-vertex(cols,:)).^2,2));
meshresolution = median(dsts); % approximation of the resolution of the mesh: median of edge lengths
clear idx rows cols dsts
nscales = round(startscale*2.^(0:2/k:(2*nscalesteps+4)/k)); % number of filtering steps
scales = meshresolution.*sqrt(2*nscales/3); % rounded values of the sigma's

%% calculate different smoothings
vertices = calculate_mesh_smoothing(vertex,faces,nscales,prior_smoothing); % calculate smoothed versions of mesh
if strcmp(type,'D')
    dextrvalue = squeeze(sqrt(sum((vertices(:,:,2:end) - vertices(:,:,1:end-1)).^2,2))); % 3D coordinates
else
    if strcmp(type,'KK')
        extrvalue = zeros(nvert,2,numel(scales));
    else
        extrvalue = zeros(nvert,numel(scales));
    end
    for i = 1:numel(scales)
        [Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertices(:,:,i),faces,options);
        switch type
            case 'K';  extrvalue(:,i) = Cgauss; % Gaussian curvature
            case 'C';  extrvalue(:,i) = Cmin.^2+Cmax.^2; % approx curvedness
            case 'H';  extrvalue(:,i) = Cmean; % Mean curvature
            case 'S';  extrvalue(:,i) = atan((Cmin+Cmax)./(Cmin-Cmax)); % approx shape index
            case 'KK';  extrvalue(:,:,i) = [Cmin Cmax]; % coordinates in principal curvature space
            otherwise
                warning('MATLAB:paramAmbiguous','type not recognised, using K (Gaussian curvature)')
                type = 'K';
                extrvalue(:,i) = Cgauss;
        end  
    end
    if strcmp(type,'KK')
        dextrvalue = squeeze(sqrt(sum((extrvalue(:,:,2:end) - extrvalue(:,:,1:end-1)).^2,2)));
    else
        dextrvalue = extrvalue(:,2:end)-extrvalue(:,1:end-1);
    end
end

maxs = zeros(nvert,2);
mins = zeros(nvert,2);
nmax = 0;
nmin = 0;
for ns = 2:numel(scales)-2;
    scale = scales(ns);
    maxsc = find(dextrvalue(:,ns)>dextrvalue(:,ns-1)&dextrvalue(:,ns)>dextrvalue(:,ns+1)); % maximum candidates
    minsc = find(dextrvalue(:,ns)<dextrvalue(:,ns-1)&dextrvalue(:,ns)<dextrvalue(:,ns+1)); % minimum candidates

    r = 1;
    for index = maxsc'
        indices = nbhmli(faces,index,r);
        indices = indices(indices~=index);
        if all(dextrvalue(index,ns)>dextrvalue(indices,ns))
            if all(dextrvalue(index,ns)>dextrvalue(indices,ns+1))
                if all(dextrvalue(index,ns)>dextrvalue(indices,ns-1))
                    nmax = nmax+1;
                    maxs(nmax,:)=[index scale];
                end
            end
        end
    end
    for index = minsc'
        indices = nbhmli(faces,index,r);
        indices = indices(indices~=index);
        if all(dextrvalue(index,ns)<dextrvalue(indices,ns))
            if all(dextrvalue(index,ns)<dextrvalue(indices,ns+1))
                if all(dextrvalue(index,ns)<dextrvalue(indices,ns-1))
                    nmin = nmin+1;
                    mins(nmin,:)=[index scale];
                end
            end
        end
    end
end
maxs = maxs(1:nmax,:);
mins = mins(1:nmin,:);
%% eliminate double keypoints: look for optimal scale
extr = [maxs;mins]; % all extrema
[un,unidx] = unique(extr(:,1)); 
unextr = extr(unidx,:); % the unique list of extrema
for index = un'
    idx = find(extr(:,1)==index);
    if numel(idx)>1 % if more than one extrema for the same vertex
        maxi = maxs(maxs(:,1)==index,:);
        mini = mins(mins(:,1)==index,:);
        for ns = 2:numel(scales)-2;
            maxi(maxi(:,2)==scales(ns),3)=ns;
            mini(mini(:,2)==scales(ns),3)=ns;
        end
        val = [2*dextrvalue(index,maxi(:,3))-dextrvalue(index,maxi(:,3)-1)-dextrvalue(index,maxi(:,3)+1) ...
              -2*dextrvalue(index,mini(:,3))+dextrvalue(index,mini(:,3)-1)+dextrvalue(index,mini(:,3)+1)]; % how good the extremum is
        [mm,i] = max(val);
        unextr(unextr(:,1)==index,2)=extr(idx(i),2); % assign the scale where the extremum is most "extreme"
    end
end
extr = unextr; % return only unique extrema
save(savepath,'extr','scales','vertices','vertex','faces')
end