

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
