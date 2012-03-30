function [n,edges] = whist(y,w,x)
% whist (weighted histogram) calculates a weighted histogram
% 
% INPUTS
% ------
% y = vector with values to make a histogram
% w = vector with weights assigned to every element of y.
% y and w must be of equal size!
% x = number of bins or an array with edges for the histogram.
%   = 8 by default
%
% OUTPUTS
% -------
% n = resulting weighted histogram
% edges = array with edges used for the histogram (identical to the input
% if edges were given
% 
% Author: Chris Maes
% 2010/04

if nargin < 3
    x = 8;
    if nargin < 2
        if nargin < 1
            error('MATLAB:hist:InvalidInput', 'You must give at least a vector y.')
        else
            if min(size(y))==1, y = y(:); end
            w = ones(size(y));
        end
    else
        if min(size(y))==1, y = y(:); end
        if min(size(w))==1, w = w(:); end
        if size(y,1)~=size(w,1)
            error('MATLAB:whist:InvalidInput', 'y and w must be of equal size')
        end
    end
end
if length(x)>1
    edges = squeeze(x);
    x = numel(edges)-1;
else
    edges = min(y):(max(y)-min(y))/x:max(y); % equally spaced bins if no edges were given
end
n = zeros(1,x);
for xi = 1:x-1
    n(xi) = sum(w(edges(xi)<=y & y<edges(xi+1))); % fill histogram
end
n(x) = sum(w(edges(x)<=y & y<=edges(x+1))); % include last boundary.