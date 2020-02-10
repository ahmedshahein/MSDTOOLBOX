function V = VectorWithoutZeros(Vector)

%
% V = VectorWithoutZeros(Vector)
%
% This function used to remove and zeros from the input vector to have
% vector with only intgers larger than zeros, which can be used as index as
% an example.
%
%   Vector: Input vector holding zeros and intgers.
%
%   Example:
%   Vector  = [0 0 3 0 5]
%   V       = [3 5]

j = 1;

for i = 1 : length(Vector),
    if Vector(i) ~= 0,
        V(j) = Vector(i);
        j=j+1;
    end
end
