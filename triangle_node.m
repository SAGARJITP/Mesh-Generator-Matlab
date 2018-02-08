%Datastructure for delaunay triangle nodes

classdef triangle_node < handle
    
    properties
        V;          %Vertex Matrix
        F;          %Opposite Face Matrix
        child_1;    %Children Node 1
        child_2;    %Children Node 1
        child_3;    %Children Node 1
        parent_1;   %Parent Node 1
        parent_2;   %Parent Node 2        
    end
    
    methods
        function triangle=triangle_node()            
            %Initiate triangle-node
            
            triangle.V=[];
            triangle.F=triangle_node.empty(3,0);           
            
            triangle.parent_1=[];
            triangle.parent_2=[];             
                      
        end
        
        %Function to create triangles based on root node and point 'p'        
        function D1=create_triangles(root,p,D)
            
            %Classify input point into triangle
            triangle=class(root,p);
            %Check if point is in interior or on edge
            [t,v]=checkonline(triangle,p);
            
            %if point is in the interior create three children triangle
            %nodes            
            if t==false              
                
                %Create and define childre node
                
                triangle.child_1=triangle_node();
                triangle.child_2=triangle_node();
                triangle.child_3=triangle_node();
                
                triangle.child_1.F(3,1)=triangle_node;
                triangle.child_2.F(3,1)=triangle_node;
                triangle.child_3.F(3,1)=triangle_node;
            
                triangle.child_1.parent_1=triangle;            
                triangle.child_2.parent_1=triangle;
                triangle.child_3.parent_1=triangle;
               
            
                triangle.child_1.V(1,:)=triangle.V(1,:);
                triangle.child_1.V(2,:)=triangle.V(2,:);
                triangle.child_1.V(3,:)=p;
                
                triangle.child_1.F(1)=triangle.child_2;
                triangle.child_1.F(2)=triangle.child_3;
                triangle.child_1.F(3)=triangle.F(3);
                                                
                sort_vertices(triangle.child_1);              
                                
                assign_opp_tri(triangle,triangle.F(3),triangle.child_1);
                
                triangle.child_2.V(1,:)=triangle.V(2,:);
                triangle.child_2.V(2,:)=triangle.V(3,:);
                triangle.child_2.V(3,:)=p;        
              
                triangle.child_2.F(1)=triangle.child_3;
                triangle.child_2.F(2)=triangle.child_1;
                triangle.child_2.F(3)=triangle.F(1);
                                
                sort_vertices(triangle.child_2);                               
               
                assign_opp_tri(triangle,triangle.F(1),triangle.child_2);
                                
                triangle.child_3.V(1,:)=triangle.V(3,:);
                triangle.child_3.V(2,:)=triangle.V(1,:);
                triangle.child_3.V(3,:)=p;        
              
                triangle.child_3.F(1)=triangle.child_1;
                triangle.child_3.F(2)=triangle.child_2;
                triangle.child_3.F(3)=triangle.F(2);
                               
                sort_vertices(triangle.child_3);                                
                
                assign_opp_tri(triangle,triangle.F(2),triangle.child_3);
                
                D(triangle==D)=[];
                D1=[D;triangle.child_1;triangle.child_2;triangle.child_3];
                                
                D1=Validedge(triangle.child_1,p,D1);
                D1=Validedge(triangle.child_2,p,D1);
                D1=Validedge(triangle.child_3,p,D1);               
                
            
            %if point is on the edge create 4 children nodes for two opposing triangles along the edge    
            elseif t==true
                
                [~,v1]=checkonline(triangle.F(v),p);   
                if v1==0
                    D1=D;
                    return;
                end
                                                             
                triangle.child_1=triangle_node();
                triangle.child_2=triangle_node();                               
                
                triangle.child_1.F(3,1)=triangle_node;
                triangle.child_2.F(3,1)=triangle_node;                                
            
                triangle.child_1.parent_1=triangle;            
                triangle.child_2.parent_1=triangle;
                
                               
                if any(any((triangle.F(v).V)))                   
                    
                    triangle.F(v).child_1=triangle_node();
                    triangle.F(v).child_2=triangle_node();
                
                    triangle.F(v).child_1.F(3,1)=triangle_node;
                    triangle.F(v).child_2.F(3,1)=triangle_node;
                
                    triangle.F(v).child_1.parent_1=triangle.F(v);
                    triangle.F(v).child_2.parent_1=triangle.F(v);
                                                                      
                end
                
                                
                triangle.child_1.V(1,:)=triangle.V(triangle_node.cw(v),:);
                triangle.child_1.V(2,:)=triangle.V(v,:);
                triangle.child_1.V(3,:)=p;        
              
                triangle.child_1.F(1)=triangle.child_2;
                triangle.child_1.F(2)=triangle.F(v).child_1;                
                triangle.child_1.F(3)=triangle.F(triangle_node.ccw(v));
                                
                sort_vertices(triangle.child_1);                             
                                                       
                assign_opp_tri(triangle,triangle.F(triangle_node.ccw(v)),triangle.child_1);
                
                triangle.child_2.V(1,:)=triangle.V(v,:);
                triangle.child_2.V(2,:)=triangle.V(triangle_node.ccw(v),:);
                triangle.child_2.V(3,:)=p;        
              
                triangle.child_2.F(1)=triangle.F(v).child_2;
                triangle.child_2.F(2)=triangle.child_1;                
                triangle.child_2.F(3)=triangle.F(triangle_node.cw(v));
              
                sort_vertices(triangle.child_2);                
                                                        
                assign_opp_tri(triangle,triangle.F(triangle_node.cw(v)),triangle.child_2);
                
                D(triangle==D)=[];
                D1=[D;triangle.child_1;triangle.child_2];                                                                      
                
                
                if any(any((triangle.F(v).V)))                                      
                
                    triangle.F(v).child_1.V(1,:)=triangle.F(v).V(v1,:);
                    triangle.F(v).child_1.V(2,:)=triangle.F(v).V(triangle_node.ccw(v1),:);
                    triangle.F(v).child_1.V(3,:)=p;
              
                    triangle.F(v).child_1.F(1)=triangle.child_1;
                    triangle.F(v).child_1.F(2)=triangle.F(v).child_2;
                    triangle.F(v).child_1.F(3)=triangle.F(v).F(triangle_node.cw(v1));
                                                
                    sort_vertices(triangle.F(v).child_1);       
                                
                
                    assign_opp_tri(triangle.F(v),triangle.F(v).F(triangle_node.cw(v1)),triangle.F(v).child_1);
                    
                    triangle.F(v).child_2.V(1,:)=triangle.F(v).V(triangle_node.cw(v1),:);
                    triangle.F(v).child_2.V(2,:)=triangle.F(v).V(v1,:);
                    triangle.F(v).child_2.V(3,:)=p;
                
                    triangle.F(v).child_2.F(1)=triangle.F(v).child_1;
                    triangle.F(v).child_2.F(2)=triangle.child_2;
                    triangle.F(v).child_2.F(3)=triangle.F(v).F(triangle_node.ccw(v1));               
                                                
                    sort_vertices(triangle.F(v).child_2);                                              
                                                
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
        
        %Sub-function to sort ertices in counter-clockwise direction
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
        
        %Sub-fucntion to assign opposite faces to triangle node
        function assign_opp_tri(triangle,opp_tri,new_tri)            
            if any(any((opp_tri.V)))
                opp_tri.F(triangle==opp_tri.F)=new_tri;                               
            end
        end
    
        %Sub-function to calculate if point is in the incirle of given
        %triangle
        function i=incircle(triangle,p)
            A=[ triangle.V(1,1),triangle.V(1,2),((triangle.V(1,1))^2 + (triangle.V(1,2))^2),1;
                triangle.V(2,1),triangle.V(2,2),((triangle.V(2,1))^2 + (triangle.V(2,2))^2),1;
                triangle.V(3,1),triangle.V(3,2),((triangle.V(3,1))^2 + (triangle.V(3,2))^2),1;
                p(1), p(2), ((p(1)^2)+(p(2))^2),1];
            
            i=det(A);
        end      
        
        %Sub-Function to flipedges if point is in the incircle of given
        %triangle
        function [t1,t2,D1]=flipedge(triangle,p,D)            
            [~,i1]=ismember(p,triangle.V,'rows');          
            
            
            [~,i2]=ismember(triangle.V(triangle_node.ccw(i1),:),triangle.F(i1).V,'rows');            
                       
            [~,v1]=checkonline(triangle.F(i1),p);
            q=triangle.F(i1).V(triangle_node.ccw(i2),:);                        
              
            triangle.child_1=triangle_node();
            triangle.child_2=triangle_node();
            
            triangle.child_1.F(3,1)=triangle_node;
            triangle.child_2.F(3,1)=triangle_node;                
                        
            triangle.child_1.parent_1=triangle;
            triangle.child_1.parent_2=triangle.F(i1);
            
            triangle.child_2.parent_1=triangle;
            triangle.child_2.parent_2=triangle.F(i1);
                       
            
            triangle.F(i1).child_1=triangle.child_1;
            triangle.F(i1).child_2=triangle.child_2;
            
            triangle.child_1.V(1,:)=p;
            triangle.child_1.V(2,:)=q;
            triangle.child_1.V(3,:)=triangle.V(triangle_node.cw(i1),:);
                                    
            triangle.child_1.F(1)=triangle.F(i1).F(i2);                           
            triangle.child_1.F(2)=triangle.F(triangle_node.ccw(i1));
            triangle.child_1.F(3)=triangle.child_2;
            
            sort_vertices(triangle.child_1);                     
            
            assign_opp_tri(triangle.F(i1),triangle.F(i1).F(i2),triangle.child_1);
            assign_opp_tri(triangle,triangle.F(triangle_node.ccw(i1)),triangle.child_1);
                                              
            triangle.child_2.V(1,:)=p;
            triangle.child_2.V(2,:)=triangle.V(triangle_node.ccw(i1),:);
            triangle.child_2.V(3,:)=q;           
                                   
            triangle.child_2.F(1)=triangle.F(i1).F(triangle_node.cw(i2));           
            triangle.child_2.F(2)=triangle.child_1;           
            triangle.child_2.F(3)=triangle.F(triangle_node.cw(i1));
            
            sort_vertices(triangle.child_2);                     
            
            assign_opp_tri(triangle.F(i1),triangle.F(i1).F(triangle_node.cw(i2)),triangle.child_2);
            assign_opp_tri(triangle,triangle.F(triangle_node.cw(i1)),triangle.child_2);
                                   
            t1= triangle.child_1;
            t2= triangle.child_2;
            
            D(triangle==D)=[];
            D(triangle.F(i1)==D)=[]; 
            D1=[D;triangle.child_1;triangle.child_2];
        end               
        
        %Sub-function to validate edge after a triangle node is created 
        function D1=Validedge(triangle,p,D)
            D1=D;
            [~,i1]=ismember(p,triangle.V,'rows');         
            
            if any(any((triangle.F(i1).V)))                
                if incircle(triangle.F(i1),p)>0                  
                    [T1,T2,D1]=flipedge(triangle,p,D);                   
                   
               
                    D1=Validedge(T1,p,D1);
                    D1=Validedge(T2,p,D1);
                end
            end
        end
        
        %Sub-function to check if point lies in the interior of given triangle 
        function t1=checkinterior(triangle,p)
            if triangle_node.side(p, triangle.V(3,:),triangle.V(1,:),triangle.V(2,:))==1 && triangle_node.side(p, triangle.V(1,:),triangle.V(2,:),triangle.V(3,:))==1 && triangle_node.side(p, triangle.V(2,:),triangle.V(3,:),triangle.V(1,:))==1
               t1=1;
            else
               t1=0;
            end           
        end
        
        %Sub-function to check if point lies on the edge of given triangle 
        function [t,v]=checkonline(triangle,p)
            t=false;
            v=0;
           
            for i=1:3
                if cross(triangle.V(triangle_node.ccw(i),:)-triangle.V(i,:), p-triangle.V(i,:))==0
                    t=true;
                    v=triangle_node.cw(i);                    
                end
            end
        end
        
        %Sub-function to classify input point to a triangle in given
        %delaunay triangulation
        function t=class(triangle,p)                   
           if checkinterior(triangle,p)==1             
              
              if isempty(triangle.child_1) && isempty(triangle.child_2) && isempty(triangle.child_3)
                  t=triangle;
                  return;
              elseif checkinterior(triangle.child_1,p)==1
                    t=class (triangle.child_1,p);
              elseif checkinterior(triangle.child_2,p)==1
                    t=class (triangle.child_2,p);
              elseif ~isempty(triangle.child_3) && checkinterior(triangle.child_3,p)==1                  
                      t=class (triangle.child_3,p);
              else
                 t=triangle;
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
    end
end


      
            
            