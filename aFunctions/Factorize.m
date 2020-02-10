function o = Factorize(x)

%
% o = Factorize(x)
%
% This function factorize the variable x, i.e. performing a long division
% procedure.
% Example:
% x = 48
% v = [2 2 2 2 3]
% 
% x = 255
% v = [3 5 17]

var =x;     % Temprory signal keeping the input value
y   = 2;    % Intial divisor
v   =[];    % Empty vector for the factors
i   = 1;    % Intial value for iteration

while x > 1,    
    if rem(x,y) == 0,
        x = x/y;
        v(i) = y;
        i = i + 1;
    else
        y = y + 1;
    end
end

if CheckForInteger(v) == 0 | v == var,
    fprintf('Non-Factorizable Number\n');
    o = nan;
else
    o = v;
end
