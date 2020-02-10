function [Fxx, f]=GetSpectrum(X, NFFT, Fs)

%
% [Fxx, f]=GetSpectrum(X, NFFT, Fs)
%
% This function to get normalized PSD of the signal X
%       X       : Signal in which we want to calculate the PSD (Power Spectral Density)
%       NFFT    : Number of points required to represent the PSD, we get it
%                 from twopowern function to satisfy that it is power of 2 points (2^n)
%       Fs      : Sampling frequency
%
% This function is imported from function developed by
% Niklas Lotze, lotze@imtek.de
%

W       = blackman(NFFT);
[ys,f]  = psd(X, NFFT, Fs, W, []);
Fxx     = ys*norm(W)^2/sum(W)^2/(Fs/NFFT);