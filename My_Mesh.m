
clear;
clc;

elem_size = 15;

%Points, faces and boundaries for the geometry
[g_points, g_faces, b] = GetInfo('circular_disc_2d.obj');


bbox = [min(g_points(:,1)), min(g_points(:,2));max(g_points(:,1)), max(g_points(:,2))];
ttol = 0.7;
delta = 0.2;

CreateMesh(ttol,bbox,elem_size,delta,g_points,b);

%% Functions

%This function creates mesh
function CreateMesh(ttol,bbox,elem_size,delta,g_points,b)

    dptol = 0.01;
    L0 = 1.2*elem_size;

    p = CreateInitialNodes(bbox,elem_size,g_points,b); %Generates initial nodes on the surface
    pfix = FixedNodes(g_points,b,elem_size);
    p = MeshPoints(p,pfix,L0);
    scatter(p(:,1),p(:,2),'x'); %Visualize initial nodes
    N = size(p,1); %No of nodes
    p_old = inf;
    
    iter = 1; %Max number of iterations
    while (iter < 500)
    
        if ( max(sqrt((p-p_old).^2),2) > ttol)
            p_old = p;
            t = delaunay(p(:,1:2)); %Delaunay triangulation
            p_mid = ((p(t(:,1),:)) + (p(t(:,2),:)) + (p(t(:,3),:)))/3;
            t = t(Distance(p_mid,g_points,b) < 0,:);
            edges = [t(:,[1,2]); t(:,[1,3]); t(:,[2,3]);]; %get all the edges of triangular mesh
            edges=unique(edges,'rows'); % Remove duplicate edges
            trimesh(t,p(:,1),p(:,2),zeros(N,1))
            view(2),axis equal,axis off,drawnow
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

%This function deletes points close to pfix and return the final p
function p = MeshPoints(p,pfix,elem_size)

    n = 1;
    for i = 1:size(pfix,1)
        for  j = 1:size(p,1)
            d = EuclideanDistance(pfix(i,:),p(j,:));
            if (d < elem_size )
               idx(n,1) = j;
               n = n+1;
            end
        end      
    end
    idx = unique(idx,'rows');
    p(idx,:)=[];
    p = [pfix;p];
end

%This function creates and returns fixed nodes on the the bounary
function pfix = FixedNodes(g_points,b,elem_size)


    pfix = g_points(b(1:size(b,1)-1),:);
    n = size(pfix,1);
    n = n+1;
    for  i = 1:(size(b,1)-1)
        d = EuclideanDistance(g_points(b(i),:),g_points(b(i+1),:));     
        if  (d/elem_size > 2)
         sub_div = floor(d/elem_size); %Get closest integer
         u_dir = Direction(g_points(b(i),:),g_points(b(i+1),:));      
            for j = 1:sub_div - 1
               pfix(n,:) = g_points(b(i),:) + j*elem_size*u_dir;
               n = n+1;
            end
        end              
    end   
end

%This function creates Initial nodes on the surface
function p = CreateInitialNodes(bbox,s,g_points,b)

   %[x,y] = meshgrid(bbox(1,1):s:bbox(2,1), bbox(1,2):s*sqrt(3/4):bbox(2,2));
   [x,y] = meshgrid(bbox(1,1):s:bbox(2,1), bbox(1,2):s:bbox(2,2));
   %x(2:2:end,:)=x(2:2:end,:)+s/2;
   p=[x(:),y(:)];
   p=[p, zeros(size(p,1),1)];
   p = p(Distance(p,g_points,b)<0,:); %keep only the points inside the boundary

end

%This function gives position of the points wrt boundary (-ve if inside the boundary)
function fd = Distance(p,g_points,b)

    fd = zeros(size(p,1),1);
    
    for i = 1:size(p,1)
        sum_pos = 0;
        sum_neg = 0;
        for j = 1:(size(b,1)-1)
            c = cross(g_points(b(j+1),:) - g_points(b(j),:), g_points(b(j),:) - p(i,:));
            if  (c(1,3) < 0)
                sum_neg = sum_neg + c;
            end
            if (c(1,3) > 0)
                sum_pos = sum_pos + c;
            end
        end
        if (sum_neg == 0 | sum_pos == 0)
            fd(i) = -1;
        else
            fd(i) = 1;
        end     
    end
end

%This function returns faces and vertices boundaries of the geometry
function [points,faces,b] = GetInfo(filename)

    Max = MaxPolygon(filename); %Get # vertices for a polygon with max vertices
    C = Polyhedron(filename,Max); %Read and store geometry information
    
    %v and f indices are stored in Index_Matrix
    Index_Matrix = C{1,1};
    
    %Numeric values from the file are stored in Data_Matrix
    Data_Matrix = C{1,2};

    [n_faces, n_vertices] = Count(Index_Matrix);
    
    points = Data_Matrix(1:n_vertices,1:3); %Get the vertices
    faces = Data_Matrix(n_vertices+1:size(Index_Matrix,1),:);
    
   
    b = boundary(points(:,1),points(:,2),0.3); %boundary vector

end

%This function gives number of vertices in polygon with highest vertices
function Max = MaxPolygon(filename)
    fileID = fopen(filename);
    D = textscan(fileID,'%s');
    fclose(fileID);
    X = D{1,1};
    N = size(X,1);

    i = 1;
    while(X{i,1} ~= 'f')
            i = i+1;
    end    

    Max = 0;
    count = 0;
    while (i <= N)
        if (X{i,1} ~= 'f')
            count = count + 1;
        end
        if (X{i,1} == 'f')
            if (count > Max)
                Max = count;
            end
            count = 0;
        end 
        i = i+1;
    end
end

%This function reads from the obj file
function C = Polyhedron(filename,Max)
    
    fmt=['%s' repmat('%f',1,Max)];
    fileID = fopen(filename);
    C = textscan( fileID,fmt,'EmptyValue',nan, 'CollectOutput', 1 );
    fclose(fileID);
    
end

%Return the no of faces and vertices for the geometry
function [n_faces,n_vertices] = Count(Index_Matrix)

    n_vertices = 0;
    i = 1;
    n = length(Index_Matrix);
    while (strcmp(Index_Matrix{i,1},'v') == 1)
            n_vertices  = n_vertices + 1;
            i=i+1;
    end
    %no_of_faces = total line - no_of_vertices
    n_faces = n - n_vertices;
end

%This function calculates Euclidean Distance between 2 points
function Eucd = EuclideanDistance(A,B)

    Eucd = sqrt( (A(1,1) - B(1,1))^2 +  (A(1,2) - B(1,2))^2);
end

%This function calculates the unit direction between two points
function m = Direction(A,B)
    m = (B - A);
    m = m/EuclideanDistance(A,B);
end

