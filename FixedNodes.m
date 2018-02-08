

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
