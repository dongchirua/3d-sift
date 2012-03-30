function [R,t,rmse] = myprocrustes(pl1,pl2);
% this function computes the Rotation (R) and translation (T) matrix from
% points in pl1 to points in pl2.
%
% ------
% INPUTS
% ------
% pl1 and pl2: structure containing
% pl1.Location= 3xN matrix with coordinates of the points from object 1
% pl2.Location= 3xN matrix with coordinates of corresponding poinst from 
%               object 2   
% 
% -------
% OUTPUTS
% -------
% R = 3x3 rotation matrix
% T = 3x1 translation matrix
% rmse = rmse for this transformation
% 
% Author: Jeroen Hermans
 
  pts1 = pl1.Location;
  pts2 = pl2.Location;

  nbPts = size(pts1,2);
  
  mean1 = mean(pts1,2);
  mean2 = mean(pts2,2);
  
  pts1C = pts1-mean1*ones(1,nbPts);
  pts2C = pts2-mean2*ones(1,nbPts);
  
  H = pts1C*pts2C';
  [U,S,V] = svd(H);
  
  R = V*U';
  
  if(det(R)<0)
    if(nbPts==3)
      V=V*diag([1 1 -1]);
      R = V*U';
    else
      error('det(R) == -1 & n > 3 : use a RANSAC alg.');
    end
  end
  
  t = mean2-R*mean1;
  
  pts1T = R*pts1 + t*ones(1,nbPts);
  err = pts2-pts1T;
  rmse = sqrt(mean(sum(err.^2)));