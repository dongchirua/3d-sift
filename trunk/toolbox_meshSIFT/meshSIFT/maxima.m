function im = maxima(x,tmax)
% this function returns the indices of all local maxima above tmax
%
% ------
% INPUTS
% ------
% x = vector
% tmax = value above which maxima must be found. (0.8*maximum by default)
%
% -------
% OUTPUTS
% -------
% im = Nx1 vector with indices of maxima
% 
% Author: Chris Maes
% 2010/04

[m,imax] = max(x);
nx = length(x);
if nargin < 2
    tmax = 0.8*m;
end
tempx = x;
i2 = imax;
ok = 0;
while ~ok
    i2 = i2+1;
    if i2 == nx+1
        i2 = 1;
        if x(i2)>x(nx)
            ok = 1;
            i2 = nx;
        end
    else
        if x(i2)>x(i2-1)
            ok = 1;
            i2 = i2-1;
        end
    end
end
i1 = imax;
ok = 0;
while ~ok
    i1 = i1-1;
    if i1 == 0
        i1 = nx;
        if x(i1)>x(1)
            ok = 1;
            i1 = 1;
        end
    else
        if x(i1)>x(i1+1)
            ok = 1;
            i1 = i1+1;
        end 
    end
end
if i2>imax
    tempx(imax:i2) = -inf;
else
    tempx(imax:end) = -inf;
    tempx(1:i2) = -inf;
end
if i1<imax
    tempx(i1:imax) = -inf;
else
    tempx(1:imax) = -inf;
    tempx(i1:end) = -inf;
end
if max(tempx) > tmax
    im = [imax maxima(tempx,tmax)];
else
    im = imax;
end