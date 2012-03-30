function [mtransf,flist1,flist2] = meshSift_posenormalisation(vertex1,faces1,vertex2,faces2)
% This function calculates the homogenous transformation matrix from mesh 1
% to mesh 2 using meshSift features.
%
% ------
% INPUTS
% ------
% vertex = Mx3 matrix,
% faces = Nx3 matrix
% 
% -------
% OUTPUTS
% -------
% mtransf = 4x4 transformation matrix from mesh 1 to mesh 2 using
%           homogenous coordinates
% flist = cell containing for every feature a structure with
%     filename = name of the file
%     kpcoordinates = coordinates of the keypoint
%     kpindex = index of the keypoint
%     kpscale = scale of the keypoint
%     orientation = canonical orientation assigned to the keypoint
%     feature = the assigned feature vector
%     numberOfFeatures = number of features assigned to this keypoint (in
%         case of multiple orientations)
% 
% Author: Chris Maes
% 2010/05

% calculate bag of features for every mesh
flist1 = meshSift(vertex1,faces1);
flist2 = meshSift(vertex2,faces2);

mtransf = calculate_MS_mtransf(flist1,flist2);

%% show meshes before posenormalisation
figure();
subs = 0; % 2 subplots in 1 figure?
if subs
    subplot(121)
end
[x,y,z] = xyz(vertex1);
trisurf(faces1,x,y,z);view(0,90);shading interp
hold on
[x,y,z] = xyz(vertex2);
trisurf(faces2,x,y,z);camlight;
axis equal;axis off
mvertex = mtransf*[vertex1 ones(size(vertex1,1),1)]';
mvertex = mvertex(1:3,:); % cut last homogenous coordinate

if subs
    subplot(122)
else
    figure;
end
[x,y,z] = xyz(mvertex');
trisurf(faces1,x,y,z);view(0,90);shading interp;
hold on
[x,y,z] = xyz(vertex2);
trisurf(faces2,x,y,z);camlight;
axis equal;axis off
end

%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%% subfunction
function [x,y,z] = xyz(n)
% extracts the three columns from an array of 3D-points
% Author: Chris Maes                                   
% 2009/09

x = n(:,1);
y = n(:,2);
z = n(:,3);
end
