function [ Q ] = bfunc(P)
M=[1 0 0 0; -3 3 0 0; 3 -6 3 0;-1 3 -3 1];

n=20; % number of values between points; 
v=linspace(0,1,n); % parameterization
u=linspace(0,1,n); % parameterization

% generate points
for i = 1:length(u)
    for j = 1:length(v)
    U = [1 u(i) u(i)^2 u(i)^3];
    V = [1 v(j) v(j)^2 v(j)^3]';
        for k = 1:3 % interpolation for each dimension of control points
            Q(i,j,k) = U * M * P(:,:,k) * M' * V;
        end
    end
end


end

