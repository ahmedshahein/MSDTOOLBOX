function e = ElementsOfVector(v)

%
% e = ElementsOfVector(v)
%
% This function used to determine the individual elements of the input vector. 
% Example:
% v = [ 1 2 5 5 5 4 8 6 2 3 7 1];
% e = [1 2 3 4 5 7 8];
% The total number of elements in v = 12 with repetitive, but e has only 7 individual elements.

h=[];
j=1;

for i = 1 : length(v)-1,
    if v(i) ~= v(i+1) & ismember(v(i),h)==0,
        h(j) = v(i);
        j=j+1;
    end
end

if ismember(v(length(v)), h),
    e = BubbleSort(h);
else
    tmp = [h v(length(v))];
    e = BubbleSort(tmp);
end
    