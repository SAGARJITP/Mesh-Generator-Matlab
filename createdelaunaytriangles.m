function T=createdelaunaytriangles(P)

P=[P(:,1:2),zeros(size(P,1),1)];

M=max(max(abs(P)));

root=triangle_node();

root.V(1,:)=[3*M,0,0];
root.V(2,:)=[0,3*M,0];
root.V(3,:)=[-3*M,-3*M,0];

root.F(3,1)=triangle_node;

root.sort_vertices;

D=[root];
k=randperm(size(P,1));

P1=P(k,:);
%P1=P;

%t=randperm(size(P1,1),1);

for r=1:size(P1,1)
    %P2=P1;
    pr=P1(r,:);
    
    D=root.create_triangles(pr,D);
end

%{%
P3=[P;root.V(1,:);root.V(2,:);root.V(3,:);0,0,0;[]];

Tri1=zeros(size(D,1),3);

for j=1:size(D,1)
    [~,l1]=ismember(D(j).V(1,:),P3, 'rows');
    [~,l2]=ismember(D(j).V(2,:),P3, 'rows');
    [~,l3]=ismember(D(j).V(3,:),P3, 'rows');
    
    Tri1(j,:)=[l1,l2,l3];
end
%}

%{
P3=[P;root.V(1,:);root.V(2,:);root.V(3,:);0,0,0;[]];
Tree=DTriTree();
%Tree.root=root;
i=1;
Tree.get_data(root,P3);


Tri1=Tree.Tri;
%}

t2=[];
n=size(P,1);

for j1=1:size(Tri1,1)
    if Tri1(j1,1)>=n+1 || Tri1(j1,2)>=n+1 || Tri1(j1,3)>=n+1
        t2=[t2;j1];
    end
end

Tri1(t2,:)=[];
%Tri1=sort(Tri1,2);
Tri2=unique(Tri1,'rows');

T=Tri2;
end
       