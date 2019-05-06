function color_surf(T,vol,p,v)

   %plots input VOL parameter over surface T triangulation with power law
   %p is how much to amplify color, between 0 and 1. Lower p's enhance
   %color more, 1 does no enhancing. v is viewing angle (optional).
   
   
   figure 

   Pts = T.Points;
   Tri = T.ConnectivityList;

   if nargin == 1
      vol = Pts(:,1);
   end

   %Normalize
   
   nvol = sign(vol).*abs(vol).^p;
   nvol = (nvol + 1)*(1/3);

   %Switch colors so red is small, blue large
   nvol = 2/3-nvol;
   patch('Faces',Tri,'Vertices',Pts,'FaceVertexCData',nvol,'FaceColor','interp','EdgeColor','none');
   daspect([1 1 1])
   axis tight

   if nargin == 4
      view(v)
   end 

   camlight('infinite')
   material dull   
   axis off

   colormap jet

end
