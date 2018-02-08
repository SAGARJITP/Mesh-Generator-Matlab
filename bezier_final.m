clear
close all
load patch_control.mat
P(:,:,:) = patch_control;

%A way of setting control points taken from user
% for k = 1:3
% P(:,:,k) = [ p00(k) p01(k) p02(k) p03(k);
%              p10(k) p11(k) p12(k) p13(k); 
%              p20(k) p21(k) p22(k) p23(k); 
%              p30(k) p31(k) p32(k) p33(k)];
% end
%conversion matrix M
M=[1 0 0 0; -3 3 0 0; 3 -6 3 0;-1 3 -3 1];

n=20; % number of values between points; 
v=linspace(0,1,n); % parameterization
u=linspace(0,1,n); % parameterization

% generate points
for i = 1:length(u)
    for j = 1:length(v)
    U = [1 u(i) u(i)^2 u(i)^3];
    V = [1 v(j) v(j)^2 v(j)^3]';
        for k = 1:3 % interpolation for each dimension of control points
            Q(i,j,k) = U * M * P(:,:,k) * M' * V;
        end
    end
end

figure 
hold on
surface(Q(:,:,1),Q(:,:,2),Q(:,:,3))
title('\bf Bezier Patch');
view(3); 

%Convert to stl b4 converting to obj
x=Q(:,:,1);
y=Q(:,:,2);
z=Q(:,:,3);
stlwrite('test.stl',x,y,z,'mode','ascii');

%set end condition to to pull stop reading from stl file
fid=fopen('test.stl','a');
fprintf(fid,'Theend\n');
fclose(fid);

%prepare for obj
fid=fopen('test.stl');
line = fgetl(fid);
i=1;
while 0==strcmp(line,'Theend')

    line = fgetl(fid);
    words = textscan(line,'%s');
    
    if strcmp(words{1,1}{1,1}, 'vertex')
      
      vx=str2double(words{1,1}{2,1});
      vy=str2double(words{1,1}{3,1});
      vz=str2double(words{1,1}{4,1});
      vertexes(i,:) = [vx, vy, vz];
      i=i+1;
    end
end
fclose(fid);

[vertex_Set,ia,ic] = unique(vertexes(:,1:3),'rows');

j=1;
for i=1:length(vertexes)
    if rem(i,3)==0
        [~,indx1]=ismember(vertexes(i-2,:),vertex_Set,'rows');
        [~,indx2]=ismember(vertexes(i-1,:),vertex_Set,'rows');
        [~,indx3]=ismember(vertexes(i,:),vertex_Set,'rows');
        faces(j,:) = [indx1, indx2, indx3];
        j=j+1;
    end
end

[rFaces,~]=size(faces);
[rVertexes,~]=size(vertex_Set);

fid2=fopen('test.obj','w');
k=1;
for i=1:(rFaces+rVertexes)
    if i <= rVertexes
        fprintf(fid2,'v ');
        fprintf(fid2,num2str(vertex_Set(i,:)));
        fprintf(fid2, '\n');
    else
        fprintf(fid2,'f ');
        fprintf(fid2,num2str(faces(k,:)));
        fprintf(fid2, '\n');
        k=k+1;
    end
    
        
end
fclose(fid2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extra stuff might need later
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if strcmp(words{1,1}{1,1}, 'vertex')
%       fprintf(fid2,'v ');
%       fprintf(fid2,words{1,1}{2,1});
%       fprintf(fid2,' ');
%       fprintf(fid2,words{1,1}{3,1});
%       fprintf(fid2,' ');
%       fprintf(fid2,words{1,1}{4,1});
%       fprintf(fid2,'\n');
%     end


