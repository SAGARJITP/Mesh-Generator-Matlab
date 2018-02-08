
%This function calculates the unit direction between two points
function m = Direction(A,B)
    m = (B - A);
    m = m/EuclideanDistance(A,B);
end