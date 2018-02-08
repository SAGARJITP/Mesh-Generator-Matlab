
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