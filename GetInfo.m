
%This function returns faces and vertices boundaries of the geometry
function [points,faces,b] = GetInfo(filename)

    Max = MaxPolygon(filename); %Get # vertices for a polygon with max vertices
    C = Polyhedron(filename,Max); %Read and store geometry information
    
    %v and f indices are stored in Index_Matrix
    Index_Matrix = C{1,1};
    
    %Numeric values from the file are stored in Data_Matrix
    Data_Matrix = C{1,2};

    [n_faces, n_vertices] = Count(Index_Matrix);
    
    points0 = Data_Matrix(1:n_vertices,1:3); %Get the vertices
    faces = Data_Matrix(n_vertices+1:size(Index_Matrix,1),:);
    
    
    
    Normalmat=zeros(n_faces,3);
    N=zeros(1,3);

    for k=1:n_faces          

        F=faces(k,:);

        v_array=zeros(3,3); %Define vector array for each vertex of given face
        %Compute vertex vectors from face definition and read vertex definition
        %from imported file:

        for k1=1:3
            p=F(1,k1);

            v_array(k1,:)=points0(p,:);           
        end

        v1=v_array(2,:)-v_array(1,:);
        v2=v_array(3,:)-v_array(2,:);
        n=computenormal(v1,v2);

        Normalmat(k,:)=n;
        N=N+n;    
    end

    Nunique=unique(Normalmat,'rows');
    N0=[0,0,1];

    if size(Nunique,1)==1 && isequal(Nunique,N0)
        points=points0;
        
    else

        Navg=N/n_faces;
        
        theta1=atan2(Navg(1),Navg(2));
        R1=[cos(theta1),-sin(theta1),0;
            sin(theta1),cos(theta1),0;
            0,0,1];

        theta2=atan2((sqrt(Navg(1)^2+Navg(2)^2)),Navg(3));
        R2=[1,0,0;
            0,cos(theta2),-sin(theta2);
            0,sin(theta2),cos(theta2)];

        %Nnew=R2*R1*Navg';

        Anew=R2*R1*points0';

        points=Anew';
    end
    
    %b = boundary(points(:,1),points(:,2),0.3); %boundary vector
    b = convhull(points(:,1),points(:,2)); %boundary vector
end
