function n=computenormal(v1,v2)
normal= [v1(2)*v2(3)-v1(3)*v2(2), v1(3)*v2(1)-v1(1)*v2(3), v1(1)*v2(2)-v1(2)*v2(1)]; %Computing normal vector using cross-product definition
normag=sqrt(dot(normal,normal));
n=normal/normag; %Unit normal vector of given vertices
end