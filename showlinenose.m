function [ output_args ] = showlinenose( xx )
%SHOWLINENOSE Summary of this function goes here
%   Detailed explanation goes here
a=imreadbw(strcat('data/bindata/',xx));
[names x, y] = textread('Nose.txt','%s %d %d');
h=find(strcmp(names,xx));
disp(h);
b=a(y(h),:);
[ymax,imax,ymin,imin] = extrema(b);
plot([1:640],b,imax,ymax,'g.',imin,ymin,'r.')
end

