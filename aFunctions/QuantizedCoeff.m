function qcoeff = QuantizedCoeff(n, Coeff_Format, coeff)

%
% qcoeff = QuantizedCoeff(n, Coeff_Format, coeff)
%
% This function used to convert fixed/floating coeff. to quantized coeff.
% in specific number of bits. It used to estimate the number of bits used
% to represent the filter coefficients, for hardware implementation. We
% tune this number of bits to maintain better IBN noise after coeff.
% quantization.
%
%       n               : number of bits for quantization
%       Coeff_Format    : define the format of input coeff. to this function, is 
%                         it fixed or floating
%       coeff           : filter coefficients to be quantized
%

q_coeff = quantizer(Coeff_Format,'floor','wrap',[n n-1]);

qcoeff  = quantize(q_coeff,coeff);