function flag = CheckForInteger(x)

%
% flag = CheckForInteger(x)
%
% This function check if the input variable or vector holds a non-integer
% value or not. If it is not integer it flags false '0' if it is integer it
% flags '1'.
% Example:
% x = 5         -> flag = 1
% x = 1.75      -> flag = 0
% x = [5 2 3]   -> flag = 1
% x = [5 2.3 3] -> flag = 0

for i = 1 : length(x),
    if rem(x(i),1) == 0,
        flag = 1;
    else
        flag = 0;
        break
    end
end