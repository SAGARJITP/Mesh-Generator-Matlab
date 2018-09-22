

%This function creates Initial nodes on the surface
function p = CreateInitialNodes(bbox,s,g_points,b)

   %[x,y] = meshgrid(bbox(1,1):s:bbox(2,1), bbox(1,2):s*sqrt(3/4):bbox(2,2));
   [x,y] = meshgrid(bbox(1,1):s:bbox(2,1), bbox(1,2):s:bbox(2,2));
    %x(2:2:end,:)=x(2:2:end,:)+s/2;
   p=[x(:),y(:)];
   p=[p, zeros(size(p,1),1)];
   p = p(Distance(p,g_points,b)<0,:); %keep only the points inside the boundary

end
