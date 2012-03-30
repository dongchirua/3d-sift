function [coords mangles match1 match2 coord1 coord2] = findcmatches(flist1,flist2,mlim)
% findmatches compares two lists of feature-structures
%
% ------
% INPUTS
% ------
% flist1 and flist2 = cell containing feature-structures
% featurelist = cell containing for every feature a structure with
%     filename = name of the file
%     kpcoordinates = coordinates of the keypoint
%     kpindex = index of the keypoint
%     kpscale = scale of the keypoint
%     orientation = canonical orientation assigned to the keypoint
%     feature = the assigned feature vector
%     numberOfFeatures = number of features assigned to this keypoint (in
%         case of multiple orientations)
% mlim = limit below which the margin must fall in order to be an
%        acceptable match (0.8 by default)
%
% -------
% OUTPUTS
% -------
% coords = Nx6 matrix containing coordinates of corresponding feature points
%          coords(i,1:3) corresponds with coords(i,4:6) (for every i)
% mangles = match angles. this vector containts for every pair in coords
%           the angle between the two featurevectors.
% match1 = vector containing for every feature in flist1 the index of the 
%          matched feature in flist2
% match2 = idem for every feature in flist1 
% margin1 = vector containing for every feature in flist1 the margin
%           between the first match in flist2 and the second match
% margin2 = idem for every feature in flist2

if nargin<3 || mlim > 1
    mlim = 0.8;
end

%% init
nf1 = numel(flist1);
nf2 = numel(flist2);
featuresize = numel(flist1{1}.feature);
mf1 = zeros(nf1,featuresize);
for i = 1:nf1
    mf1(i,:) = flist1{i}.feature;
end
mf2 = zeros(nf2,featuresize);
for i = 1:nf2
    mf2(i,:) = flist2{i}.feature;
end

%% calculations
% angles = zeros(nf2,nf1);
% for it_dim=1:featuresize
%     angles = angles + (ones(nf2,1)*mf1(:,it_dim)'  - mf2(:,it_dim)*ones(1,nf1) ).^2;
% end
% angles = sqrt(angles);
angles = real(acos((repmat(1./sqrt(sum(mf2.^2,2)),1,featuresize).*mf2)*(repmat(1./sqrt(sum(mf1.^2,2)),1,featuresize).*mf1)'));

%% matches for flist1
[min1,idx] = min(angles,[],1);
inds = (0:size(angles,2)-1)*size(angles,1)+idx;
anglestemp = angles;
anglestemp(inds)=inf;
min2 = min(anglestemp,[],1);
margin1 = min1./min2;
match1 = zeros(1,nf1);
match1(margin1<mlim) = idx(margin1<mlim);

coord1 = zeros(nf1,6);
mangle1 = zeros(1,nf1);
nc = 0;
for i = 1:nf1
    if match1(i)~=0
        mcoord = [flist1{i}.kpcoordinates flist2{match1(i)}.kpcoordinates];
        ind = strfind(reshape(coord1',1,[]),mcoord);
        if isempty(ind(rem(ind,6) == 1))
            nc = nc+1;
            coord1(nc,:) = mcoord;
            mangle1(nc) = min1(i);
        end
    end
end
coord1 = coord1(1:nc,:);
mangle1 = mangle1(1:nc);

%% matches for flist2
anglestemp = angles.';
[min1,idx2] = min(anglestemp,[],1);
inds = (0:size(angles,1)-1)*size(angles,2)+idx2;
anglestemp(inds)=inf;
min2 = min(anglestemp,[],1);
margin2 = min1./min2;
match2 = zeros(1,nf2);
match2(margin2<mlim) = idx2(margin2<mlim);

coord2 = zeros(nf2,6);
mangle2 = zeros(1,nf2);
nc = 0;
for i = 1:nf2
    if match2(i)~=0
        mcoord = [flist2{i}.kpcoordinates flist1{match2(i)}.kpcoordinates];
        ind = strfind(reshape(coord2',1,[]),mcoord);
        if isempty(ind(rem(ind,6) == 1))
            nc = nc+1;
            coord2(nc,:) = mcoord;
            mangle2(nc) = min1(i);
        end
    end
end
coord2 = coord2(1:nc,:);
mangle2 = mangle2(1:nc);

coord = [coord1;[coord2(:,4:6) coord2(:,1:3)]];
coords = zeros(size(coord));
mangles = zeros(1,size(coord,1));
nc = 0;
for i = 1:size(coord,1)
    ind = strfind(reshape(coord(i+1:end,:)',1,[]),coord(i,:));
    if ~isempty(ind(rem(ind,6) == 1))
        nc = nc+1;
        coords(nc,:) = coord(i,:);
        if i<=size(coord1,1)
            mangles(nc) = mangle1(i);
        else
            mangles(nc) = mangle2(i-size(coord1,1));
        end
    end
end
coords = coords(1:nc,:);
mangles = mangles(1:nc);
end