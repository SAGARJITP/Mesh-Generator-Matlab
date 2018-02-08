

%This function reads from the obj file
function C = Polyhedron(filename,Max)
    
    fmt=['%s' repmat('%f',1,Max)];
    fileID = fopen(filename);
    C = textscan( fileID,fmt,'EmptyValue',nan, 'CollectOutput', 1 );
    fclose(fileID);
    
end
