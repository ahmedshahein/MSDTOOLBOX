function z = FindVector(x, y)

%
% z = FindVector(x, y)
%
% This function allocate the index from vectr x if vector y belongs to x.
%   Example:
%   x=[1 3 5 7 9 11 13];
%   y=[5 7];
%   z = find_vector(x, y);
%   z = 3 4

width = length(y)-1;
k = 1;%
z = [];
for i = 1 : length(x)-width,
    if x(i:i+width) == y,
        z(k,:) = i:i+width;
        k = k + 1;
    end
end