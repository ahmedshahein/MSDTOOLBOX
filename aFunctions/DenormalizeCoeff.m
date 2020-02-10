function dcoeff = DenormalizeCoeff(coeff)

% This function denormalize the normalized flter coefficients, to have
% filter response normalized to zero, i.e. the response at dc will be at
% 0dB.

dcoeff = coeff./sum(coeff);