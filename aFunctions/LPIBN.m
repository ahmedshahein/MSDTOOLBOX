function [IBN,Sig] = LPIBN(data_vector, Fsignal, Fs, Fb, plot_psd)

%
% [IBN,Sig] = LPIBN(data_vector, Fsignal, Fs, Fb, plot_psd)
%
% This function used to setimate the IBN for the input vector of the
% function
%
%   data_vector: Vector holds the output pattern of the DSM or the Decimator
%   fsignal: Signal frequency
%   Fs: Sampling frequency
%   OSR: Oversampling frequency
%   plot_psd: Optinal flag to plot the PSD of the function input signal or not, 1 -> plot 0 -> don't plot
%
% This funcion based on the function developed by Maurits Ortmanns
% 'Maurits.Ortmanns@imtek.de' called 'FileSpektrumAuto_simu.m'
%

% Turn off the warnign of log of zero
log10(0);
[msg, msgid] = lastwarn;
s = warning('off', msgid);

fsig=Fsignal;

amp=data_vector;

N=length(amp);
power=0;
while N>2^power,
   power=power+1;
end;
Np=2^(power-1);

W=blackman(Np);
[ys,f]=psd(amp,Np,Fs,W,[]);
sp_a=ys*norm(W)^2/sum(W)^2*4;
sp_psd=ys*norm(W)^2/sum(W)^2/(Fs/Np);
sp_psd_tmp=sp_psd;

% Signal peak bin and signal cutoff bin
fsignr=round(fsig/(Fs/(Np)));
fstoppnr=round(Fb/(Fs/(Np)));

% Signal width bins
n=1;
vergl=sp_psd_tmp(fsignr);
while (min(sp_psd_tmp((fsignr+n):(fsignr+n+2)))<vergl),
   vergl=sp_psd_tmp(fsignr+n);
   n=n+1;
   if n==fsignr,
      n=n-2;
      break;
   end;
end;

% Temporary signal without signal
sp_ersatz_tmp = (sp_psd_tmp(fsignr-(n))+sp_psd_tmp(fsignr+(n)))/2;

% Removing the signal from the pattern
for i=0:n-1,
    sp_psd_tmp(fsignr-i)=sp_ersatz_tmp;
    sp_psd_tmp(fsignr+i)=sp_ersatz_tmp;
end;

% The IBN is computed now the achievement of the individual peaks
sp_psd_bin=sp_psd_tmp*(Fs/Np);              % Noise power density
IBN_tmp=0;
%startnr=fsignr;
startnr=1; % The strat point for calculating the IBN
IBN_tmp=sum(sp_psd_bin(startnr:(fstoppnr+1)));

%The signal and the IBN in dB
Sig=10*log10(sum(sp_psd(fsignr-n:fsignr+n)*(Fs/Np)));
IBN=10*log10(IBN_tmp);

%Plot the PSD with and without the signal
if (nargin == 5 & plot_psd == 1),
    figure
    semilogx(f,10*log10(sp_psd),'k-',f,10*log10(sp_psd_tmp));
    grid on;
    axis auto
    xlabel('Frequeny- Hz');
    ylabel('PSD - dB')
    title(['PSD IBN=', num2str(IBN),', Sig=', num2str(Sig)]);
end





