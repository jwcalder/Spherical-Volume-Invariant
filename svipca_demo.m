%Demo script for computing principal curvatures 
%via PCA on local neighborhoods and plotting 
%on Stanford dragon
close all;
clearvars;
addpath('c_code');
load('meshes/dragon.mat');

%TR = T; %High resolution
TR = T_d1; %Medium resolution
%TR = T_d2; %Low resolution

[S,K1,K2] = svipca(TR,1);

%Plotting
p = 0.6;

%K1
color_surf(TR,K1,p,[-15 15]);
title('K1');

%K2
color_surf(TR,K2,p,[-15 15]);
title('K2');

%Gauss curvature
color_surf(TR,K1.*K2,p,[-15 15]);
title('Gauss Curvature');
