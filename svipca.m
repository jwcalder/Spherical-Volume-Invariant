function [Sout,K1,K2,V1,V2,V3] = svipca(TR,R,ID)
%  Computes spherical volume invariant using method from
%
%  "Computation of circular area and spherical volume invariants via boundary integrals",
%  Riley O'Neill, Pedro Angulo-Umana, Jeff Calder, Bo Hessburg, Peter Olver, Cheri Shakiban, and Katrina Yezzi-Woodley, preprint, 2019.
%
%  Usage 
%
%     [Sout,K1,K2,V1,V2,V3] = svipca(TR,R,ID)
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
%     K1 = Largest principal curvature
%     K2 = Smallest principal curvature
%     V1 = Direction associated with K1
%     V2 = Direction associated with K2
%     V3 = Surface normal (orthogonal to span(V1,V2))
%
%  Uses the mex c code svipca_mex.c. Use mexmake.m to compile
%
%  Authors: Riley O'Neill, Jeff Calder, 2019

   prog = true;   %Toggle progress bar on or off (code is faster with prog=false)
   eps_svi = 1;     %Error tolerance for refinement integration (set eps=100 or larger for no refinement)
   eps_pca = 1;     %Error tolerance for refinement integration (set eps=100 or larger for no refinement)

   P = TR.Points;
   T = TR.ConnectivityList;
   len = length(P);
   if nargin == 3
       if (any(ID>1) || length(ID)~=len)   %if input is point indices form
          id = zeros(length(TR.Points),1);
          id(ID) = 1;
          ID = logical(id);
      end
   else 
      ID = logical(ones(length(P),1));
   end

   Sout = zeros(length(R),length(P));
   Mout = zeros(length(R),length(P),3,3);
   for i=1:length(R)
      r = R(i);
      [S,M] = svipca_mex(P',int32(T'-1),r,logical(ID),eps_svi,eps_pca,prog);
      M = reshape(M',[length(P),3,3]);
      Sout(i,:) = S;
      Mout(i,:,:,:) = M;
   end
   [K1,K2,V1,V2,V3] = MtoK(TR,Mout,Sout,R,logical(ID)); 
   Sout = Sout';
end

function [K1,K2,V1,V2,V3] = MtoK(T,M,vsr,r,ID)
% triangulation T, PCA matrices M and vsr (vol) directly from svipca, radii r. 
% converts PCA matrices to principal curvatures K1&K2, 
% principal directions V1, V2, & V3.


   nn = length(ID);
   rlen = length(r);
   Norm = vertexNormal(T);

   %initialization:
   K1 = zeros(nn,rlen);
   K2 = K1;
   V1 = zeros(nn,3*rlen);
   V2 = V1;
   V3 = V2;

   II = [1:nn]'; %original pt indices
   I = II(ID); %points we want to define outputs for

   for k=1:rlen
      
      lim2 = 3*k;     %lim1 & lim2 ranges for output vectors
      lim1 = lim2 - 2;
      l = lim1:lim2;
      
      L1 = zeros(nn,1);
      L2 = zeros(nn,1);
      L3 = zeros(nn,1);
      
      for j=1:length(I)
          i = I(j);
         A = squeeze(M(k,i,:,:)); 
         [V,D] = eig(A);
         D = diag(D); 

         %Check which eigenvector is normal to surface:
         a = zeros(3,1);      
         a(1) = Norm(i,:)*V(:,1); 
         a(2) = Norm(i,:)*V(:,2);
         a(3) = Norm(i,:)*V(:,3);
         [m j] = max(abs(a));
         %so all outward
       
         switch j
            case 1
               L1(i) = D(2);
               L2(i) = D(3);
               L3(i) = D(1);
               
               V1(i,l) = V(:,2);
               V2(i,l) = V(:,3);
               V3(i,l) = V(:,1);
            case 2
               L1(i) = D(1);
               L2(i) = D(3);
               L3(i) = D(2);
               
               V1(i,l) = V(:,1);
               V2(i,l) = V(:,3);
               V3(i,l) = V(:,2);
            case 3
               L1(i) = D(1);
               L2(i) = D(2);
               L3(i) = D(3);
               
               V1(i,l) = V(:,1);
               V2(i,l) = V(:,2);
               V3(i,l) = V(:,3);
         end
      end
      
      %using eig's only
      Kdiff = (L1-L2)*24/(pi*r(k)^6);
      Ksum = (16*pi*r(k)^3/3 - 8*squeeze(vsr(k,:)))/(pi*r(k)^4);
      Ksum = Ksum';
      k1t = (Kdiff + Ksum)./2;
      k2t = (Ksum - Kdiff)./2;
      
      %want to ensure k1>k2:
      J = k1t > k2t; %logical
      K1(:,k) = J.*k1t + (1-J).*k2t; %if k1 max, keep it as k1, else swap
      K2(:,k) = (1-J).*k1t + J.*k2t; 
      v1t = V1(:,l); 
      v2t = V2(:,l);
      V1(:,l) = J.*v1t + (1-J).*v2t; %so V1 corresponds to K1
      V2(:,l) = (1-J).*v1t + J.*v2t;
       
      %now for quality control: if volume is not defined:
      visnegative = vsr(k,:) == -1;
      K1(visnegative,k) = 0; 
      K2(visnegative,k) = 0;
      V1(visnegative,l) = 0;
      V2(visnegative,l) = 0;
      V3(visnegative,l) = 0;
      
      vecneg = sum(V3(:,l).*Norm<0,2)<0;
      V3(vecneg,l) = -V3(vecneg,l);
      V2(vecneg,l) = -V2(vecneg,l);
      V1(vecneg,l) = -V1(vecneg,l);
      
      %implementing right hand rule:
      rhr = sum(V3(:,l).*cross(V1(:,l),V2(:,l)),2) < 0;
      V1(rhr,l) = -V1(rhr,l);

   end 
end 
%function [K1,K2,V1,V2,V3] = MtoK(T,M,vsr,r)
% triangulation T, PCA matrices M and vsr (vol) directly from svipca, radii r. 
% converts PCA matrices to principal curvatures K1&K2, (K1>K2)
% principal directions V1, V2, & V3.
%
%
%   n = length(T.Points);
%
%   Norm = vertexNormal(T);
%
%   for k=1:length(r)
%      lim2 = 3*k;     %lim1 & lim2 ranges for output vectors
%      lim1 = lim2 - 2;
%      l = lim1:lim2;
%      
%      L1 = zeros(n,1);
%      L2 = zeros(n,1);
%      L3 = zeros(n,1);
%      
%      for i=1:n
%         A = squeeze(M(k,i,:,:)); 
%         [V,D] = eig(A);
%         D = diag(D); 
%
%         
%         %Check which eigenvector is normal to surface:
%         a = zeros(3,1);
%
%         
%         a(1) = Norm(i,:)*V(:,1); 
%         a(2) = Norm(i,:)*V(:,2);
%         a(3) = Norm(i,:)*V(:,3);
%         [m j] = max(abs(a));
%         
%         %so all outward
%        
%    
%         switch j
%            case 1
%               L1(i) = D(2);
%               L2(i) = D(3);
%               L3(i) = D(1);
%               
%               V1(i,l) = V(:,2);
%               V2(i,l) = V(:,3);
%               V3(i,l) = V(:,1);
%            case 2
%               L1(i) = D(1);
%               L2(i) = D(3);
%               L3(i) = D(2);
%               
%               V1(i,l) = V(:,1);
%               V2(i,l) = V(:,3);
%               V3(i,l) = V(:,2);
%            case 3
%               L1(i) = D(1);
%               L2(i) = D(2);
%               L3(i) = D(3);
%               
%               V1(i,l) = V(:,1);
%               V2(i,l) = V(:,2);
%               V3(i,l) = V(:,3);
%         end
%      end
%      
%      %using eig's only
%      Kdiff = (L1-L2)*24/(pi*r(k)^6);
%      Ksum = (16*pi*r(k)^3/3 - 8*squeeze(vsr(k,:)))/(pi*r(k)^4);
%      Ksum = Ksum';
%      k1t = (Kdiff + Ksum)./2;
%      k2t = (Ksum - Kdiff)./2;
%      
%      %want to ensure k1>k2:
%      J = k1t > k2t; %logical
%      K1(:,k) = J.*k1t + (1-J).*k2t; %if k1 max, keep it as k1, else swap
%      K2(:,k) = (1-J).*k1t + J.*k2t; 
%      v1t = V1(:,l); 
%      v2t = V2(:,l);
%      V1(:,l) = J.*v1t + (1-J).*v2t; %so V1 corresponds to K1
%      V2(:,l) = (1-J).*v1t + J.*v2t;
%       
%      %now for quality control: if volume is not defined:
%      visnegative = vsr(k,:) == -1;
%      K1(visnegative,k) = 0; 
%      K2(visnegative,k) = 0;
%      V1(visnegative,l) = 0;
%      V2(visnegative,l) = 0;
%      V3(visnegative,l) = 0;
%      
%      vecneg = sum(V3(:,l).*Norm<0,2)<0;
%      V3(vecneg,l) = -V3(vecneg,l);
%      V2(vecneg,l) = -V2(vecneg,l);
%      V1(vecneg,l) = -V1(vecneg,l);
%      
%      %implementing right hand rule:
%      rhr = sum(V3(:,l).*cross(V1(:,l),V2(:,l)),2) < 0;
%      V1(rhr,l) = -V1(rhr,l);
%       
%   end 
%end 
%function [K1,K2,V1,V2,V3] = MtoK(T,M,vsr,r,ID)
%   %triangulation T, M and vsr (vol) directly from svipca, radii r. 
%
%   %K1<K2
%
%   n = length(T.Points);
%
%   Norm = vertexNormal(T);
%   K1 = zeros(n,length(r));
%   K2 = zeros(n,length(r));
%   V1 = zeros(n,3*length(r));
%   V2 = zeros(n,3*length(r));
%   V3 = zeros(n,3*length(r));
%
%   for k=1:length(r)
%      lim2 = 3*k;     %lim1 & lim2 ranges for output vectors
%      lim1 = lim2 - 2;
%      
%      L1 = zeros(n,1);
%      L2 = zeros(n,1);
%      L3 = zeros(n,1);
%      
%      for i=1:n
%         if ID(i)
%            A = squeeze(M(k,i,:,:)); 
%            [V,D] = eig(A);
%            D = diag(D); %eigenvalues
%
%            
%            %Check which eigenvector is normal to surface:
%            a = zeros(3,1);
%            a(1) = abs(Norm(i,:)*V(:,1)); 
%            a(2) = abs(Norm(i,:)*V(:,2));
%            a(3) = abs(Norm(i,:)*V(:,3));
%            [m j] = max(a);
%            
%            switch j
%               case 1
%                  L1(i) = D(2);
%                  L2(i) = D(3);
%                  L3(i) = D(1);
%                  
%                  V1(i,lim1:lim2) = V(:,2);
%                  V2(i,lim1:lim2) = V(:,3);
%                  V3(i,lim1:lim2) = V(:,1);
%               case 2
%                  L1(i) = D(1);
%                  L2(i) = D(3);
%                  L3(i) = D(2);
%                  
%                  V1(i,lim1:lim2) = V(:,1);
%                  V2(i,lim1:lim2) = V(:,3);
%                  V3(i,lim1:lim2) = V(:,2);
%               case 3
%                  L1(i) = D(1);
%                  L2(i) = D(2);
%                  L3(i) = D(3);
%                  
%                  V1(i,lim1:lim2) = V(:,1);
%                  V2(i,lim1:lim2) = V(:,2);
%                  V3(i,lim1:lim2) = V(:,3);
%            end
%         end
%      end
%      
%      %using eig's only
%      Kdiff = (L1-L2)*24/(pi*r(k)^6);
%      Ksum = (16*pi*r(k)^3/3 - 8*squeeze(vsr(k,:)))/(pi*r(k)^4);
%      Ksum = Ksum';
%      k1t = (Kdiff + Ksum)./2;
%      k2t = (Ksum - Kdiff)./2;
%      
%      %want to ensure k1>k2:
%      J = k1t > k2t; %logical
%      K1(:,k) = J.*k1t + (1-J).*k2t; %if k1 max, keep it as k1, else swap
%      K2(:,k) = (1-J).*k1t + J.*k2t; 
%      v1t = V1(:,lim1:lim2); 
%      v2t = V2(:,lim1:lim2);
%      V1(:,lim1:lim2) = J.*v1t + (1-J).*v2t; %so V1 corresponds to K1
%      V2(:,lim1:lim2) = (1-J).*v1t + J.*v2t;
%       
%      %now for quality control: if volume is not defined:
%      isnegative = vsr(k,:) == -1;
%      K1(isnegative,k) = 0; 
%      K2(isnegative,k) = 0;
%      V1(isnegative,lim1:lim2) = 0;
%      V2(isnegative,lim1:lim2) = 0;
%      V3(isnegative,lim1:lim2) = 0;
%      
%     
%   end 
%end 

