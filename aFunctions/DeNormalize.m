function coeffUn = DeNormalize(coeff, Q)

% 
% coeffUn = DeNormalize(coeff, Q)
%
% This function denormalize the normalized filter coefficients.
%
%   coeff:  Normalized filter coefficients
%   Q:      Quantization bit width, it has to be the same Q used by the
%           function 'NormalizedCoeff'.
%
%   coeffUn:Denormalized filter coefficients.
%
% P.S.
% As a check up you can sum(coeffUn), it has to be equal or less than 1.

b = coeff./(2^(Q-1)-1);
coeffUn = b./sum(b);