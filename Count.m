

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