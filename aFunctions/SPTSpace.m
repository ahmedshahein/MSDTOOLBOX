function NSPT_Vector = SPTSpace(Nb, n)

% NSPT_Vector = SPTSpace(Nb, n)
% This function exports the dense SPT space.
% 
%   Nb:     Maximum number of non-zero digitas for SPT factor
%   n:      Maximum power of two factor
%
% Example:
% n = 5
% if Nb = 2 then
% 1SPT = 1 2 4 8 16 
% 2SPT = 3 6 12 24
% so the SPT space = 1 2 3 4 6 8 12 16 24
%
% if Nb = 3 then
% 1SPT = 1 2 4 8 16 
% 2SPT = 3 6 12 24
% 3SPT = 7 14 28
% so the SPT space = 1 2 3 4 6 7 8 12 14 16 24 28

NSPT_Matrix = zeros(Nb, n);

for i = 1 : Nb,
    lengths(i) = n-(i-1);
end

for j = 1 : Nb,
    NSPT_Matrix(j,1:lengths(j)) = GenerateSPTs(j, n);
end
    
k=1;
for i = 1 : n,
    for j = 1 : Nb,
        columns(k) = i-j+1;
        if i >= Nb,
            rows(k) = j;
        elseif i == 1 & j == 1,
            rows(k) = 1;
        elseif i ~= 1 & i < Nb & j <= i,
            rows(k) = j;
        else
            rows(k) = 0;
        end
        k = k+1;
    end
end

col = [];
row = [];
m = 1;
n = 1;
for l = 1 : prod(size(NSPT_Matrix)),
    if columns(l) <= 0,
        m = m;
    else
        col(m) = columns(l);
        m = m + 1;
    end
    if rows(l) <= 0,
        n = n;
    else
        row(n) = rows(l);
        n = n + 1;
    end
end

for o = 1 : length(row),
    NSPT_Vector(o) = NSPT_Matrix(row(o), col(o));
end

%%%%%%%%%%%%%%%
% SUBFUNCTION %
%%%%%%%%%%%%%%%
function NPT = GenerateSPTs(Nb, n)

% This function export the Nb-SPT factor numbers for Nb up to n bits.
%
%   Nb:     Maximum number of non-zero digitas for SPT factor
%   n:      Maximum power of two factor
%
% Example:
% n = 10 -> maximum power-of-two factor is 2^10-1 = 512
% Nb = 1 -> each coefficient only has 1 power-of-two factor
%           1 2 4 8 16 == 2^0 2^1 2^2 2^3 2^4
% Nb = 2 -> each coefficient has 2 power-of-two factors
%           3 6 12 24  == (1+2) (2+4) (4+8) (8+16)
% Nb = 3 -> each coefficient has 3 power-of-two factors
%           7 14 28    == (1+2+4) (2+4+8) (4+8+16) 

if Nb > n,
    fprintf('WARNING: It is recommended to make n larger than or equal Nb');
end

SPT = [];
for i = 1 : n,
    SPT(i) = 2^(i-1);
end

NPT = [];
for j = 1 : n-Nb+1,
    NPT(j) = sum(SPT(j:j+Nb-1));
end