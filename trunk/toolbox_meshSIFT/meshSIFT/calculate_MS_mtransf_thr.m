function [mtransf,rmse] = calculate_MS_mtransf_thr(flist1,flist2,thr)
% this function calculates the homogenous transformation matrix from mesh 1
% to mesh 2, based upon the list of features from every mesh.
%
% ------
% INPUTS
% ------
% flist1 and flist2 = cell containing feature-structures of mesh 1 and 2
% featurelist = cell containing for every feature a structure with
%     filename = name of the file
%     kpcoordinates = coordinates of the keypoint
%     kpindex = index of the keypoint
%     kpscale = scale of the keypoint
%     orientation = canonical orientation assigned to the keypoint
%     feature = the assigned feature vector
%     numberOfFeatures = number of features assigned to this keypoint (in
%         case of multiple orientations)
%     thr = threshold
%
% -------
% OUTPUTS
% -------
% mtransf = 4x4 transformation matrix from mesh 1 to mesh 2 using
%           homogenous coordinates
% find corresponding scale-space-extrema
if nargin<3
    thr = 0.8;
end
stop = 0;
coords = [];
while size(coords,1) < 5 && ~stop
    coords = findcmatches(flist1,flist2,thr);
    % calculate transformation matrix
    pl1.Location = coords(:,1:3).';
    pl2.Location = coords(:,4:6).';
    if size(coords,1) < 5
        thr = thr + (1-thr)/2;
        if thr==1
            stop = 1;
        end
        warning(['Only ' int2str(size(coords,1)) ' correspondence(s) found, threshold increased to ' num2str(thr)])
    end
end
disp([int2str(size(coords,1)) ' correspondences found'])
params.QuantilePercentage = 0.5;
params.NbTrialPoints = 4;
params.NbTrials = 100;

try
    % RANSAC transformation matrix
    [R,t,rmse] = RANSAC_LeastQuantileOfSquareErrors(pl1,pl2,params);
    mtransf = [R t;0 0 0 1];
catch
    mtransf = eye(4); rmse = Inf;
    %error 'Ransac failed'
end
end