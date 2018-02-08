
%This function calculates Euclidean Distance between 2 points
function Eucd = EuclideanDistance(A,B)

    Eucd = sqrt( (A(1,1) - B(1,1))^2 +  (A(1,2) - B(1,2))^2);
end
