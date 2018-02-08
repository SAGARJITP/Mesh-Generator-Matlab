

%This function creates mesh
function CreateMesh(ttol,bbox,elem_size,delta,g_points,b,obj)


    dptol = 0.01;
    L0 = 1.2*elem_size;

    p = CreateInitialNodes(bbox,elem_size,g_points,b); %Generates initial nodes on the surface
    pfix = FixedNodes(g_points,b,elem_size);
    p = MeshPoints(p,pfix,elem_size);
    axes(obj);
    scatter(p(:,1),p(:,2),'x'); %Visualize initial nodes
    N = size(p,1); %No of nodes
    p_old = inf;
    
    iter = 1; %Max number of iterations
    while (iter < 100)
    
        if ( max(sqrt((p-p_old).^2),2) > ttol)
            p_old = p;
            %t = createdelaunaytriangles(p(:,1:3)); %Delaunay triangulation
            t = delaunay(p(:,1:2)); %Delaunay triangulation
            p_mid = ((p(t(:,1),:)) + (p(t(:,2),:)) + (p(t(:,3),:)))/3;
            t = t(Distance(p_mid,g_points,b) < 0,:);
            edges = [t(:,[1,2]); t(:,[1,3]); t(:,[2,3]);]; %get all the edges of triangular mesh
            edges=unique(edges,'rows'); % Remove duplicate edges
            axes(obj);
            trimesh(t,p(:,1),p(:,2),zeros(N,1))
            %triplot(handles.Mesh,t,p(:,1),p(:,2))
            view(2),axis equal,axis on,drawnow
        end 
        
        edge_vec = (p(edges(:,1),:) - p(edges(:,2),:)); %edge_vec stors the vector for each edge
        L = sqrt(sum(edge_vec.^2,2)); %calculate length of each edges
        F = (L0-L);
        %F = max(L0-L,0); %Calculte the magnitude of force on each bar
        %F(1:size(pfix,1),:)=0;
        F_vec = F.*edge_vec./(L);
        
        
        %update the position of p
        for i = 1:size(edge_vec,1)  
            if (edges(i,1) > size(pfix,1))
                 p(edges(i,1),:) =  p(edges(i,1),:) + delta*F_vec(i,:);  
            end
                  
        end    
        
        if max(sqrt(sum(delta*F_vec.^2,2)))<0.035
            break;
        end
        
        iter = iter + 1;
    end
    
    
end