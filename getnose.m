function [x,y ] = getnose( t )
%GETNOSE Summary of this function goes here
%   Detailed explanation goes here
z=t;
zz=z(180:360,120:size(z,1)-120);
zmax=max(max(zz));
[r,c] = find(zz==zmax);
r=r+180;
c=c+120;
x=r(1)
y=c(1)
end

