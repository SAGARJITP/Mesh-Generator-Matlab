

%This function calculates Euclidean Distance between 2 points in an array
function Eucd = EuclideanDistance_arr(A,B)
    Eucd = sum((A-B).^2).^0.5;
end