function [x,y,z] = xyz(n)
% extracts the three columns from an array of 3D-points
% Author: Chris Maes                                   
% 2009/09

x = n(:,1);
y = n(:,2);
z = n(:,3);
end