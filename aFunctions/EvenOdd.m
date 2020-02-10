function flag = EvenOdd(x)

%
% flag = EvenOdd(x)
%
% This function to verify whether the number is and Even (flag = 1) or Odd
% (flag = 0).

if mod(x, 2) == 0,
    flag = 1;
else
    flag = 0;
end