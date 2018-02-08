
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