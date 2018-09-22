classdef triangle_node < handle
    
    properties
        V;       %Vertex 1
        F;     
        child_1;
        child_2;
        child_3;
        parent_1;
        parent_2;
        leaf;
    end
    
    methods
        function triangle=triangle_node()%(p1,p2,p3)
            %triangle.V=zeros(3,3);
            triangle.V=[];
            triangle.F=triangle_node.empty(3,0);           
            
            triangle.parent_1=[];
            triangle.parent_2=[];  
            triangle.leaf=0;
                      
        end
                                
        function D1=create_triangles(root,p,D)
            triangle=class(root,p);
            [t,v]=checkonline(triangle,p);
            
            if t==false              
                
                triangle.child_1=triangle_node();
                triangle.child_2=triangle_node();
                triangle.child_3=triangle_node();
                
                triangle.child_1.F(3,1)=triangle_node;
                triangle.child_2.F(3,1)=triangle_node;
                triangle.child_3.F(3,1)=triangle_node;
            
                triangle.child_1.parent_1=triangle;            
                triangle.child_2.parent_1=triangle;
                triangle.child_3.parent_1=triangle;
                
                triangle.leaf=1;
            
                triangle.child_1.V(1,:)=triangle.V(1,:);
                triangle.child_1.V(2,:)=triangle.V(2,:);
                triangle.child_1.V(3,:)=p;        
              
                %triangle.child_1.F=[triangle.child_2;triangle.child_3;triangle.F(3);];
                
                triangle.child_1.F(1)=triangle.child_2;
                triangle.child_1.F(2)=triangle.child_3;
                triangle.child_1.F(3)=triangle.F(3);
                                                
                sort_vertices(triangle.child_1);              
                
                %assign_opp_tri(triangle.child_1,p);
                assign_opp_tri(triangle,triangle.F(3),triangle.child_1);
                
                triangle.child_2.V(1,:)=triangle.V(2,:);
                triangle.child_2.V(2,:)=triangle.V(3,:);
                triangle.child_2.V(3,:)=p;        
              
                triangle.child_2.F(1)=triangle.child_3;
                triangle.child_2.F(2)=triangle.child_1;
                triangle.child_2.F(3)=triangle.F(1);
                                
                sort_vertices(triangle.child_2);                
               
                %assign_opp_tri(triangle.child_2,p);
                assign_opp_tri(triangle,triangle.F(1),triangle.child_2);
                                
                triangle.child_3.V(1,:)=triangle.V(3,:);
                triangle.child_3.V(2,:)=triangle.V(1,:);
                triangle.child_3.V(3,:)=p;        
              
                triangle.child_3.F(1)=triangle.child_1;
                triangle.child_3.F(2)=triangle.child_2;
                triangle.child_3.F(3)=triangle.F(2);
                               
                sort_vertices(triangle.child_3);
                                
                %assign_opp_tri(triangle.child_3,p);
                assign_opp_tri(triangle,triangle.F(2),triangle.child_3);
                
                D(triangle==D)=[];
                D1=[D;triangle.child_1;triangle.child_2;triangle.child_3];
                                
                D1=Validedge(triangle.child_1,p,D1);
                D1=Validedge(triangle.child_2,p,D1);
                D1=Validedge(triangle.child_3,p,D1);               
                
            
            elseif t==true
                %if v~=0
                    
                %adj_triangle=triangle.F(v);
                                             
                triangle.child_1=triangle_node();
                triangle.child_2=triangle_node();                               
                
                triangle.child_1.F(3,1)=triangle_node;
                triangle.child_2.F(3,1)=triangle_node;                                
            
                triangle.child_1.parent_1=triangle;            
                triangle.child_2.parent_1=triangle;
                
                triangle.leaf=1;                            
                             
                
                if any(any((triangle.F(v).V)))                   
                    
                    triangle.F(v).child_1=triangle_node();
                    triangle.F(v).child_2=triangle_node();
                
                    triangle.F(v).child_1.F(3,1)=triangle_node;
                    triangle.F(v).child_2.F(3,1)=triangle_node;
                
                    triangle.F(v).child_1.parent_1=triangle.F(v);
                    triangle.F(v).child_2.parent_1=triangle.F(v);
                
                    triangle.F(v).leaf=1;                                    
                end
                
                                
                triangle.child_1.V(1,:)=triangle.V(triangle_node.cw(v),:);
                triangle.child_1.V(2,:)=triangle.V(v,:);
                triangle.child_1.V(3,:)=p;        
              
                triangle.child_1.F(1)=triangle.child_2;
                triangle.child_1.F(2)=triangle.F(v).child_1;                
                triangle.child_1.F(3)=triangle.F(triangle_node.ccw(v));
                                
                sort_vertices(triangle.child_1);                             
                
                %assign_opp_tri(triangle.child_1,p);                              
                assign_opp_tri(triangle,triangle.F(triangle_node.ccw(v)),triangle.child_1);
                
                triangle.child_2.V(1,:)=triangle.V(v,:);
                triangle.child_2.V(2,:)=triangle.V(triangle_node.ccw(v),:);
                triangle.child_2.V(3,:)=p;        
              
                triangle.child_2.F(1)=triangle.F(v).child_2;
                triangle.child_2.F(2)=triangle.child_1;                
                triangle.child_2.F(3)=triangle.F(triangle_node.cw(v));
              
                sort_vertices(triangle.child_2);                
                
                %assign_opp_tri(triangle.child_2,p);                         
                assign_opp_tri(triangle,triangle.F(triangle_node.cw(v)),triangle.child_2);
                
                D(triangle==D)=[];
                D1=[D;triangle.child_1;triangle.child_2];                           
                                               
                %{
                [~,i]=ismember(triangle.V(triangle_node.ccw(v),:),triangle.F(v).V,'rows');
                p3=triangle.F(v).V(triangle_node.ccw(i));
                %}
                
                if any(any((triangle.F(v).V)))
                                        
                    
                    [~,v1]=checkonline(triangle.F(v),p);                    
                    %i=triangle_node.ccw(v1);
                    %[~,i]=ismember(triangle.V(triangle_node.ccw(v),:),triangle.F(v).V,'rows');
                    %v1=triangle_node.ccw(i);
                    
                    %v=zeros(3,1);
                    %{
                    for i=1:3
                        [~,z]=checkonline(triangle.F(i),p);                       
                        v(i)=z;
                    end
                    %}
                    %v1=max(v);
                    %i1=triangle_node.cw(v1);
                    %q=triangle.V(i1,:);
                    
                    %q=triangle.F(v).V((triangle==triangle.F(v).F),:);
                    %[~,v1]=ismember(q,triangle.F(v).V,'rows');
                                                              
                %adj_triangle.child_1.V(1,:)=triangle.V(triangle_node.cw(v),:);
                %{
                triangle.F(v).child_1.V(1,:)=triangle.F(v).V(triangle_node.cw(i),:);
                triangle.F(v).child_1.V(2,:)=p3;
                triangle.F(v).child_1.V(3,:)=p;        
                %}
                
                    triangle.F(v).child_1.V(1,:)=triangle.F(v).V(v1,:);
                    triangle.F(v).child_1.V(2,:)=triangle.F(v).V(triangle_node.ccw(v1),:);
                    triangle.F(v).child_1.V(3,:)=p;
              
                    triangle.F(v).child_1.F(1)=triangle.child_1;
                    triangle.F(v).child_1.F(2)=triangle.F(v).child_2;
                    triangle.F(v).child_1.F(3)=triangle.F(v).F(triangle_node.cw(v1));
                                                
                    sort_vertices(triangle.F(v).child_1);            
                                
                %assign_opp_tri(triangle.F(v).child_1,p);
                    assign_opp_tri(triangle.F(v),triangle.F(v).F(triangle_node.cw(v1)),triangle.F(v).child_1);
                                                
                %adj_triangle.child_2.V(1,:)=triangle.V(triangle_node.ccw(v),:);
                %{
                triangle.F(v).child_2.V(1,:)=triangle.F(v).V(i,:);
                triangle.F(v).child_2.V(2,:)=p3;
                triangle.F(v).child_2.V(3,:)=p;        
                %}
              
                    triangle.F(v).child_2.V(1,:)=triangle.F(v).V(triangle_node.cw(v1),:);
                    triangle.F(v).child_2.V(2,:)=triangle.F(v).V(v1,:);
                    triangle.F(v).child_2.V(3,:)=p;
                
                    triangle.F(v).child_2.F(1)=triangle.F(v).child_1;
                    triangle.F(v).child_2.F(2)=triangle.child_2;
                    triangle.F(v).child_2.F(3)=triangle.F(v).F(triangle_node.ccw(v1));               
                                                
                    sort_vertices(triangle.F(v).child_2);
                                              
                %assign_opp_tri(triangle.F(v).child_2,p);                                
                    assign_opp_tri(triangle.F(v),triangle.F(v).F(triangle_node.ccw(v1)),triangle.F(v).child_2);
                
                    D1(triangle.F(v)==D1)=[];
                    D1=[D1;triangle.F(v).child_1;triangle.F(v).child_2];
                                               
                    D1=Validedge(triangle.F(v).child_1,p,D1);
                    D1=Validedge(triangle.F(v).child_2,p,D1);                    
                end
                
                D1=Validedge(triangle.child_1,p,D1);
                D1=Validedge(triangle.child_2,p,D1);               
                %end
                
            end          
        end
        
        function sort_vertices(triangle)
            
            M1=[triangle.V(1,1),triangle.V(1,2),1;
                triangle.V(2,1),triangle.V(2,2),1;
                triangle.V(3,1),triangle.V(3,2),1];
            
            if det(M1)<0
                t=triangle.V(2,:);
                triangle.V(2,:)=triangle.V(3,:);
                triangle.V(3,:)=t;
                
                f=triangle.F(2);
                triangle.F(2)=triangle.F(3);
                triangle.F(3)=f;
            end
        end
        
        %{
        function assign_opp_tri(triangle,p)            
            [~,i1]=ismember(p,triangle.V,'rows');
            
            %if ~isempty(triangle.F(i1).V)
            if any(any((triangle.F(i1).V)))
                %opp_triangle=triangle.F(i1);
                [~,i2]=ismember(triangle.V(triangle_node.ccw(i1),:),triangle.F(i1).V,'rows');
                triangle.F(i1).F(triangle_node.ccw(i2))=triangle;
            end
        end
        %}
        
        function assign_opp_tri(triangle,opp_tri,new_tri)            
            if any(any((opp_tri.V)))
                %i=find(opp_tri.F==triangle);
                %{
                i1=0;
                for i=1:3
                    if isequal(triangle,opp_tri.F(i))
                        opp_tri.F(i)=new_tri;
                        i1=i;
                    end
                end
                %}
                %z=find((triangle==opp_tri.F)==1);
                %i=triangle_node.cw(z);
                
                opp_tri.F(triangle==opp_tri.F)=new_tri;
                
                %i=triangle_node.cw(i1);
            end
        end
    
        function i=incircle(triangle,p)
            A=[ triangle.V(1,1),triangle.V(1,2),((triangle.V(1,1))^2 + (triangle.V(1,2))^2),1;
                triangle.V(2,1),triangle.V(2,2),((triangle.V(2,1))^2 + (triangle.V(2,2))^2),1;
                triangle.V(3,1),triangle.V(3,2),((triangle.V(3,1))^2 + (triangle.V(3,2))^2),1;
                p(1), p(2), ((p(1)^2)+(p(2))^2),1];
            
            i=det(A);
        end
        
        %{
        function [p0,p1,p2,tri]=pointcheck(triangle,p)
            if triangle.vertex_1==p
                p0=triangle.vertex_1;
                p1=triangle.vertex_2;
                p2=triangle.vertex_3;
                tri=triangle.opp_tri_1;
                
            elseif triangle.vertex_2==p
                p0=triangle.vertex_2;
                p1=triangle.vertex_3;
                p2=triangle.vertex_1;
                tri=triangle.opp_tri_2;
                
            elseif triangle.vertex_3==p
                p0=triangle.vertex_3;
                p1=triangle.vertex_1;
                p2=triangle.vertex_2;
                tri=triangle.opp_tri_3;
            end
        end      
        %}
        
        %{
        function q=facecheck(triangle,adj_triangle,p)
            [~,~,p2,~]=pointcheck(triangle,p);
                        
            if adj_triangle.vertex_1==p2               
                q=adj_triangle.vertex_3;
                                
            elseif adj_triangle.vertex_2==p2
                q=adj_triangle.vertex_1;
                               
            elseif adj_triangle.vertex_3==p2                               
                q=adj_triangle.vertex_2;
            end
            
            %{
            if triangle1.opp_tri_1==triangle2
                p=triangle1.vertex_1;
            elseif triangle1.opp_tri_2==triangle2
                p=triangle1.vertex_2;
            elseif triangle1.opp_tri_3==triangle2
                p=triangle1.vertex_3;
            end
            %}
        end
          %}       
        
        function [t1,t2,D1]=flipedge(triangle,p,D)            
            [~,i1]=ismember(p,triangle.V,'rows');
            
            %{
            if any(any((triangle.F(i1).V)))
                i2=0;
                %i=find(opp_tri.F==triangle);
                for i=1:3                    
                    if isequal(triangle,triangle.F(i1).F(i))
                        i2=triangle_node.cw(i);
                        q=triangle.F(i1).V(i,:);
                    end
                end
            end
            %}
         
            
            %triangle.F(i1)
            %x=triangle.F(i1);
            
            %[~,i2]=ismember(triangle.V(triangle_node.ccw(i1),:),triangle.F(i1).V,'rows');
            
            %if i2~=0            
            %[~,v1]=checkonline(triangle.F(i1),p);
            %q=triangle.F(i1).V(triangle_node.ccw(i2),:);
            
            q=triangle.F(i1).V((triangle==triangle.F(i1).F),:);
            [~,i]=ismember(q,triangle.F(i1).V,'rows');
            i2=triangle_node.cw(i);
                       
            triangle.child_1=triangle_node();
            triangle.child_2=triangle_node();
            
            triangle.child_1.F(3,1)=triangle_node;
            triangle.child_2.F(3,1)=triangle_node;                
                        
            triangle.child_1.parent_1=triangle;
            triangle.child_1.parent_2=triangle.F(i1);
            
            triangle.child_2.parent_1=triangle;
            triangle.child_2.parent_2=triangle.F(i1);
            
            triangle.leaf=1;
            triangle.F(i1).leaf=1;
            
            triangle.F(i1).child_1=triangle.child_1;
            triangle.F(i1).child_2=triangle.child_2;
            
            triangle.child_1.V(1,:)=p;
            triangle.child_1.V(2,:)=q;
            triangle.child_1.V(3,:)=triangle.V(triangle_node.cw(i1),:);
                                    
            triangle.child_1.F(1)=triangle.F(i1).F(i2);                           
            triangle.child_1.F(2)=triangle.F(triangle_node.ccw(i1));
            triangle.child_1.F(3)=triangle.child_2;
            
            sort_vertices(triangle.child_1);
                        
            %assign_opp_tri(triangle.child_1,p);       
            %assign_opp_tri(triangle.child_1,q);            
            
            assign_opp_tri(triangle.F(i1),triangle.F(i1).F(i2),triangle.child_1);
            assign_opp_tri(triangle,triangle.F(triangle_node.ccw(i1)),triangle.child_1);
                                              
            triangle.child_2.V(1,:)=p;
            triangle.child_2.V(2,:)=triangle.V(triangle_node.ccw(i1),:);
            triangle.child_2.V(3,:)=q;           
                                   
            triangle.child_2.F(1)=triangle.F(i1).F(triangle_node.cw(i2));           
            triangle.child_2.F(2)=triangle.child_1;           
            triangle.child_2.F(3)=triangle.F(triangle_node.cw(i1));
            
            sort_vertices(triangle.child_2);                        
            
            %assign_opp_tri(triangle.child_2,p);                         
            %assign_opp_tri(triangle.child_2,q);
            
            assign_opp_tri(triangle.F(i1),triangle.F(i1).F(triangle_node.cw(i2)),triangle.child_2);
            assign_opp_tri(triangle,triangle.F(triangle_node.cw(i1)),triangle.child_2);
                                   
            t1= triangle.child_1;
            t2= triangle.child_2;
            
            D(triangle==D)=[];
            D(triangle.F(i1)==D)=[]; 
            D1=[D;triangle.child_1;triangle.child_2];
            %{
            else
                t1= triangle;
                t2= triangle.F(i1);                
            end
            %}
        end
        
        %{
        function tri=find_opp_tri(triangle,p)
            if triangle.vertex_1==p               
                tri=triangle.opp_tri_1;
                
            elseif triangle.vertex_2==p                
                tri=triangle.opp_tri_2;
                
            elseif triangle.vertex_3==p                
                tri=triangle.opp_tri_3;
            end
        end
           %}
                        
        function D1=Validedge(triangle,p,D)
            D1=D;
            [~,i1]=ismember(p,triangle.V,'rows');
            %adj_tri=triangle.F(i1);
                        
            %if ~isempty(triangle.F(i1).V)
            if any(any((triangle.F(i1).V)))                
                if incircle(triangle.F(i1),p)>0                  
                    [T1,T2,D1]=flipedge(triangle,p,D);
               
                    D1=Validedge(T1,p,D1);
                    D1=Validedge(T2,p,D1);                    
                end
            end
        end
                        
        function t1=checkinterior(triangle,p)            
            if triangle_node.side(p, triangle.V(3,:),triangle.V(1,:),triangle.V(2,:))==1 && triangle_node.side(p, triangle.V(1,:),triangle.V(2,:),triangle.V(3,:))==1 && triangle_node.side(p, triangle.V(2,:),triangle.V(3,:),triangle.V(1,:))==1
               t1=1;
            else
               t1=0;
            end           
        end
        
        function [t,v]=checkonline(triangle,p)
            t=false;
            v=0;
           
            for i=1:3
                if cross(triangle.V(triangle_node.ccw(i),:)-triangle.V(i,:), p-triangle.V(i,:))==0
                    t=true;
                    v=triangle_node.cw(i);
                    %break;
                end
            end
        end
        
        function t=class(triangle,p)                   
           if checkinterior(triangle,p)==1
              
              %{
              exist triangle.child_1; A1=ans; %,'triangle_node');
              exist triangle.child_2; A2=ans; %'triangle_node');
              exist triangle.child_3; A3=ans; %,'triangle_node');
              %}
              
              if isempty(triangle.child_1) && isempty(triangle.child_2) && isempty(triangle.child_3)
                  t=triangle;
                  return;
              elseif checkinterior(triangle.child_1,p)==1
                    t=class (triangle.child_1,p);
              elseif checkinterior(triangle.child_2,p)==1
                    t=class (triangle.child_2,p);
              elseif ~isempty(triangle.child_3) && checkinterior(triangle.child_3,p)==1
                  %if checkinterior(triangle.child_3,p)==1
                      t=class (triangle.child_3,p);              
                  %end
              %else
              %   t=triangle;
              end            
          end
        end
             
    end
    
    methods (Static)
        
        function t2=side(p1,p2,a,b)
            cp1= cross(b-a,p1-a);
            cp2= cross(b-a,p2-a);            
            
            if dot(cp1,cp2)>=0
                  t2=1;
            else
                  t2=0;
            end            
        end
        
        function i2=ccw(i)            
            if i~=0
                i1=mod(i,3);
                i2=i1+1;
            end
        end
        
        function i2=cw(i)
            if i~=0
                i1=mod(i-2,3);
                i2=i1+1;
            end
        end
        %}
    end
end


      
            
            