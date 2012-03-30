% example_meshSIFT.m
% ------------------
%
% This file contains a simple example to find and match features between
% two scans of the same person.

clear all
close all
clc

addpath(genpath('../meshSIFT/'))
addpath(genpath('../toolbox_fast_marching/'))
addpath(genpath('../toolbox_graph/'))
%% loading
s1 = load('Dirk_cr_ds_1');
s2 = load('Dirk_cr_ds_2');

%% compute features
fl1 = meshSift(s1.vertex,s1.faces);
fl2 = meshSift(s2.vertex,s2.faces);

%% match features
coords = findcmatches(fl1,fl2);

%% visualization
offset = 15;
figure;
trisurf(s1.faces,s1.vertex(:,1),s1.vertex(:,2),s1.vertex(:,3),ones(length(s1.vertex),1)); 
hold on
trisurf(s2.faces,s2.vertex(:,1)+offset,s2.vertex(:,2),s2.vertex(:,3),ones(length(s2.vertex),1)); 
for k=1:size(coords,1)
    plot3([coords(k,1) coords(k,4)+offset],[coords(k,2) coords(k,5)],[coords(k,3) coords(k,6)],'LineWidth',1.5);
end
plot3(coords(:,1),coords(:,2),coords(:,3),'og','MarkerSize',4,'MarkerFaceColor','g')
plot3(coords(:,4)+offset,coords(:,5),coords(:,6),'or','MarkerSize',4,'MarkerFaceColor','r')
axis equal; shading interp; lighting phong; camlight; axis off; view(0,90)
colormap gray