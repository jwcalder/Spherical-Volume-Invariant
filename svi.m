function [Sout, Sgam] = svi(TR,R,ID)
%  Computes spherical volume invariant using method from
%
%  "Computation of circular area and spherical volume invariants via boundary integrals",
%  Riley O'Neill, Pedro Angulo-Umana, Jeff Calder, Bo Hessburg, Peter Olver, Cheri Shakiban, and Katrina Yezzi-Woodley, preprint, 2019.
%
%  Usage 
%
%     [Sout,Sgam] = svi(TR,R,ID);
%  
%  where
%
%     TR = surface triangulation
%     R = list of radii to compute invariant
%     ID = optional boolean vector of vertices or list of vertex indices 
%          at which to compute volume
%
%  Output
%     
%     Sout(i,j) = Spherical volume invariant at vertex i for radius j
%     Sgam = Gamma (see paper)
%
%  Uses the mex c code svi_mex.c. Use mexmake.m to compile
%
%  Authors: Riley O'Neill, Jeff Calder, 2019

   prog = true;   %Toggle progress bar on or off (code is faster with prog=false)
   eps = 1;    %Error tolerance for refinement integration (set eps=100 or larger for no refinement)

   P = TR.Points;
   T = TR.ConnectivityList;
   len = length(P);

   if nargin == 3
      if (any(ID>1) || length(ID)~=len)   %if input is point indices
          id = zeros(length(TR.Points),1);
          id(ID) = 1;
          ID = id;
      end
   else
      ID = ones(length(P),1);
   end

   Sout = zeros(length(R),length(P));
   for i=1:length(R)
      r = R(i);
      [S,Sgam] = svi_mex(P',int32(T'-1),r,logical(ID),eps,prog);
      Sout(i,:) = S;
   end
  
   Sout = Sout';
   Sgam = Sgam';
end
