function B = replace(A,n1,n2)
% looks in a matrix A for occurrences of number n1 and replaces them with
% n2.
% Author: Chris Maes
% 2009/10

[nrows,ncols] = size(A);
B = A;
for i = 1: nrows
    for j = 1:ncols
        if A(i,j) == n1
            B(i,j) = n2;
        end
    end    
end
end