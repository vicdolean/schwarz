function [M,S] = iteration_matrix(n,a,b)

d1 = b*ones(2*n-3,1);
d1(2:2:end) = 0;
d2 = a*ones(2*n-4,1); 
d2(2:2:end) = 0;
d3 = a*ones(2*n-4,1)-d2;

M = diag(d1,1)+diag(d1,-1)+diag(d2,2)+diag(d3,-2);
S = eig(M);