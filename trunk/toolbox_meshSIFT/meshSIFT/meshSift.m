function featurelist = meshSift(vertex,faces,filename)
% calculates scale space extrema and calculates features for those extrema
% uses: ssextrema.m and calculate_ssfeatures.m
%
% INPUTS
% ------
% vertex = Mx3 matrix,
% face = Nx3 matrix
% filename (optional) = name assigned to every feature
%
% OUTPUT
% ------
% featurelist = cell containing for every feature a structure with
% filename = name of the file
% kpname (optional) = name of the keypoint
% kpcoordinates = coordinates of the keypoint
% kpindex = index of the keypoint
% kpscale = scale of the keypoint
% orientation = canonical orientation assigned to the keypoint
% feature = the assigned feature vector
% numberOfFeatures = number of features assigned to this keypoint (in
%                    case of multiple orientations)
%
% Author: Chris Maes
% 2010/04

%% init
if nargin < 3
    filename = 'unknown';
end

%% calculate scale space extrema

%%% default parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% extrparams.k = 4;
% extrparams.type = 'H';
% extrparams.options.curvature_smoothing = 3;
% extrparams.nscalesteps = 5;
% extrparams.startscale = 5; 1?
% extrparams.prior_smoothing = 1;
% [extr,scales,vertices] = ssextrema(vertex,faces,extrparams);
[extr,scales,vertices] = ssextrema(vertex,faces);

%% calculate features of the scale space extrema
% parameters
featparams.scales = scales;
featparams.vertices = vertices;
featparams.filename = filename;

%%% default parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% featparams.patchmode = 1; % manner of distributing patches around the point of intrest: 1 = square, 2 = logpolar
% featparams.orientationmode = 1; % mode to calculate the normal (based on normals or on Umax)
% featparams.fscale = 3; % factor to scale overall size of the feature
% featparams.forient = 3; % factor for the size of the neighborhood concerned in defining the orientation.
% featparams.fsigorient = 1.5; % factor for the sigma of the weighting factors for defining the orientation. (1/1.5?)
% featparams.fsighist = 1.5; % factor for the sigma of the weighting factors for the histogram (1/0.5?)
% featparams.fsighistangle = 0.5; % factor for the sigma for the angle of the weighting factors for the histogram
% featparams.nbinsorient = 360; % # bins for the histogram to select the orientation.
% featparams.nbinsfeature = 8; % # bins for the histogram in the feature vector.
% featparams.fdist = 1.5; % factor for the distance from the keypoint, where to center the local patch for the histogram.
% featparams.fpatch = 1.25; % factor for the radius for the local patch of the histogram
% featparams.normmode = 1; % mode to normalise the histograms (mode 2 means renormalising...)
% featparams.curvature_smoothing = 3; % size of smoothing ring for curvature calculation

featurelist=calculate_meshSiftFeatures(vertex,faces,extr,featparams);