%Demo script for computing spherical volume invariant 
%and plotting on the Stanforddragon
close all;
clearvars;
addpath('c_code');
load('meshes/dragon.mat');

%TR = T; %High resolution
TR = T_d1; %Medium resolution
%TR = T_d2; %Low resolution

%Compute spherical volume invariant
S = svi(TR,1); 

%Plot
color_surf(TR,S,1,[-15,15]);
title('Spherical Volume Invariant');
