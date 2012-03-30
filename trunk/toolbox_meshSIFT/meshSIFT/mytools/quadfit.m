function coef = quadfit(list)
% This function generates a quadratic surface from coefficients in coef and
% evaluates these for the datapoints in data.
% coef = [a b c d e f]
% data = [x y z]
% function: ax^2+by^2+cxy+dx+ey+f = z
% Author: Chris Maes
% 2009/09

[nrows,ncols] = size(list);
if ncols ~= 3
    warning('MATLAB:paramAmbiguous','3D-points are needed as input')
    coef = 0;
    return
elseif nrows < 9
    warning('MATLAB:paramAmbiguous','At least 9 points are needed to define a unique quadratic surface')
    coef = 0;
    return
end
A = ones(nrows,6);
A(:,1:2) = list(:,1:2).^2; % x^2 and y^2
A(:,3) = list(:,1).*list(:,2); % xy
A(:,4) = list(:,1); % x
A(:,5) = list(:,2); % y
b = list(:,3); % z
coef = inv(A.'*A)*A.'*b;
end