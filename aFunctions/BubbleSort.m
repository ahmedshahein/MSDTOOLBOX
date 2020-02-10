function x = BubbleSort(x)

%
% x = BubbleSort(x)
%
% This function rearrange the input vector in ascending order.
% Example:
% x = [4 7 2 5 9];
% x = BubbleSort(x)
%   = [2 4 5 7 9]

for i =1:length(x)-1,
    for j = 1:length(x)-1,
        if x(j) > x(j+1),
            temp=x(j);
            x(j)=x(j+1);
            x(j+1)=temp;
        end
    end
end