function [ME MSE] = FilterCoefficientMeanError(qCoeff, nCoeff, Fs)

%
% [ME MSE] = FilterCoefficientMeanError(qCoeff, nCoeff)
%
% This function calculates the mean error (ME) and the mean square error
% between the quantized coefficients (qCoeff) and the scaled or normalized
% coefficients (nCoeff).

[Hq Fq]=freqz(qCoeff./sum(qCoeff), 1, 2^10, Fs);
Magq = 20*log10(abs(Hq));


[Hn Fn]=freqz(nCoeff./sum(nCoeff), 1, 2^10, Fs);
Magn = 20*log10(abs(Hn));

figure
plot(Fn, Magn)
hold on
plot(Fq, Magq, 'r')

MSE = mean( abs(Hn - Hq).^2 );
ME = mean( abs(Hn - Hq) );