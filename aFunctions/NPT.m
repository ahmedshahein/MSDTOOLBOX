function out = NPT(in, Nb, n)

% out = NPT(in, N, n)
% This function exports the correspondence Signed-Power-of-Two (SPT) value for
% the input value 'in'.
%
%   in :    Input integer
%   Nb:      Number of ones in for the desired SPT, N belongs to (1, 2, 3) only
%   n:      Quantization bit width for this number

out2=[];
single_PT = SPTSpace(2,n);
double_PT = SPTSpace(Nb,n);
  
    
    for l = 1 : length(in),
        for m = 1 : length(double_PT)-1,
            if in(l) > 0,
                if in(l) >= double_PT(m) & in(l) <= double_PT(m+1),
                    diff1 = double_PT(m+1) - in(l);
                    diff2 = in(l) - double_PT(m);
                    if min([diff1 diff2]) == diff1,
                        out2(l) = double_PT(m+1);
                    else
                        out2(l) = double_PT(m);
                    end
                end
            else
                if abs(in(l)) >= double_PT(m) & abs(in(l)) <= double_PT(m+1),
                    diff1 = double_PT(m+1) - abs(in(l));
                    diff2 = abs(in(l)) - double_PT(m);
                    if min([diff1 diff2]) == diff1,
                        out2(l) = -1*double_PT(m+1);
                    else
                        out2(l) = -1*double_PT(m);
                    end
                elseif in(l) == 0,
                    out2(l) = 0;
                end   
            end
        end
    end 

    if Nb == 1,
        out = SPT(in);
    else,
        out = out2;
    end