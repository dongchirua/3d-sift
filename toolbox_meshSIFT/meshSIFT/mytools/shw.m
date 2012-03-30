function shw(n,holdon,varargin)
% quickly shows a scatterplot of an array of 3D-points
% Author: Chris Maes                                  
% 2009/09

if nargin < 3
    style = '.';
    if nargin < 2
        holdon = 0;
    end
end

if holdon == 0
    figure()
else
    hold on
end
scatter3(n(:,1),n(:,2),n(:,3),varargin{:})
axis equal
hold off
end