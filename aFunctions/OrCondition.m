function z = OrCondition(a,b,c)

%
% z = OrCondition(a,b,c)
%
% This function to test if the number 'c' belongs to the numbers from 'a'
% to 'b', if yes the output is '1' if not the output will be '0'.
%
%       a : left/low number range
%       b : right/high number range
%       c : number to check if it belongs to the specified number range or not
%       z : output, 1/0 yes/no
%

% Definting the nymber range
m = (a:b);

% Examining if the inout 'c' belongs to this number range or not
for i = 1:length(m),
    if m(i)==c,
        x(i) = 1;
    else
        x(i) = 0;
    end
end

x;

if sum(x)==1,
    z=1;
else
    z=0;
end