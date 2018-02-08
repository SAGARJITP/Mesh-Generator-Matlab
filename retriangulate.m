function P1=retriangulate(filename,P,P0)

%{
A=importdata('surface.obj');
%}

%P(b,:)=[];

A=importdata(filename);

l=length(A.textdata); %Determine length of data file
vertices=0; 
faces=0;

% Compute no. of vertices iteratively by checking for character 'v' in each
% line of imported file

for j1=1:l
    if (strcmp(A.textdata(j1,1),'v'))
        vertices=vertices+1;
    end
end

% Compute no. of vertices iteratively by checking for character 'f' in each
% line of imported file

for j2=1:l
    if (strcmp(A.textdata(j2,1),'f'))
        faces=faces+1;
    end
end

Facedata(faces,1)=triangle_node;
Facedata1(faces,1)=triangle_node;

for k=1:faces          
                   
    F=A.data(vertices+k,:);
               
    for k1=1:3
        p0=F(1,k1);
               
        Facedata(k).V(k1,:)=A.data(p0,:);
        Facedata1(k).V(k1,:)=P0(p0,:); 
    end
end

P1=zeros(size(P,1),3);

t=triangle_node;
t1=triangle_node;

for j=1:size(P,1)
    p=[P(j,1:2),0];
    for j1=1:faces        
        
        if checkinterior1([Facedata1(j1).V(:,1:2),[0;0;0]],p)==1
            t1=Facedata1(j1);       
            t=Facedata(j1);
        end
    end
    
    
    a=(t.V(2,:)-t.V(1,:))/norm(t.V(2,:)-t.V(1,:));
    b=(t.V(3,:)-t.V(1,:))/norm(t.V(3,:)-t.V(1,:));
    
    %{%
    a1=(t1.V(2,:)-t1.V(1,:))/norm(t1.V(2,:)-t1.V(1,:));
    b1=(t1.V(3,:)-t1.V(1,:))/norm(t1.V(3,:)-t1.V(1,:));
    %}
    
    %{
    a1=(t1.V(2,1:2)-t1.V(1,1:2))/norm(t1.V(2,1:2)-t1.V(1,1:2));
    b1=(t1.V(3,1:2)-t1.V(1,1:2))/norm(t1.V(3,1:2)-t1.V(1,1:2));
    %}
    
    %p1=p(:,1:2)-t1.V(1,1:2);
    p1=p-t1.V(1,:);
    
    C=[a1(1),b1(1);
        a1(2),b1(2)];
    
    T=inv(C)*[p1(1);p1(2)];
    u=T(1);
    v=T(2);
    %u=dot(p(1,1:2),a1(1,1:2));
    %v=dot(p(1,1:2),b1(1,1:2));
    
    P1(j,:)=(u*a+v*b)+t.V(1,:);    
end

%{
Aorig=A.data(1:vertices,:)';

figure;
scatter3(P1(:,1),P1(:,2),P1(:,3));

figure;
scatter3(Aorig(1,:),Aorig(2,:),Aorig(3,:));
%}
end

    
