function [deviation_in_coeff norm_coeff mixed_coeff count] = OptiMixedCoeff(coeff, Q, Nb, deviation_value)

%
% [deviation_in_coeff norm_coeff mixed_coeff count] = OptiMixedCoeff(coeff, Q, deviation_value)
%
% This function exports the mixed coefficients vector due to acceptable
% deviation value. The mixed vector holds coefficients rounded to SPT and
% NPT and normal normalized coefficients. How it works:
% If the deviation of the coeff. is less than or equal the deviation value,
% the coeff. will be rounded to the value having the smallest deviation,
% i.e.
% coeff. = 53 -> SPT = 64, NPT = 48
% 64 - 53 = 11, while 53 - 48 = 5
% so the coeff. will be rounded to 48 not 64
%

count = 0;

norm_coeff = NormalizedCoeff(coeff, Q);
spt_coeff = SPT(norm_coeff);
npt_coeff = NPT(norm_coeff, Nb, Q);

dev_spt = abs(norm_coeff - spt_coeff);
dev_npt = abs(norm_coeff - npt_coeff);

for i = 1 : length(coeff),
    tmp_dev     = [dev_spt(i) dev_npt(i)];
    tmp_coeff   = [spt_coeff(i) npt_coeff(i)];
    if (tmp_dev(1) <= deviation_value & tmp_dev(2) <= deviation_value),
        if min(tmp_dev),
            mixed_coeff(i) = tmp_coeff(IndexOfMinValue(tmp_coeff));
            deviation_in_coeff(i) = min(tmp_dev);
            count = count + 1;
        else
            mixed_coeff(i) = norm_coeff(i);
            deviation_in_coeff(i) = 0;
        end
    else
        mixed_coeff(i) = norm_coeff(i);
        deviation_in_coeff(i) = 0;
    end
end