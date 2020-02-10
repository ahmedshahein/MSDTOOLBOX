function RT = ComputationalEffort(Fs, OSR, delta_F, K, M, rp, rc, print_RT)

%
% RT = ComputationalEffort(Fs, OSR, delta_F, K, M, rp, rc, print_RT)
% 
% This function used to estimate the computational effort of FIR filter taking 
% into consideration the effect of; number of decimation stages, decimation 
% factor for each stage overall decimation factor, sampling frequency,
% passband and cutoff band ripples.
%
%   Fs:                    Sampling frequency 
%   OSR:                   Over Sampling Ration
%   delta_F:               Decimator transtion bandwidth
%   K:                     Number of decimation stages
%   M:                     Decimation vector for decimation factor at each stage
%   rp:                    Pass-band ripples
%   rc:                    Cutoff-band ripples
%   print_RT:              Flag to print RT on the command window
%
%   RT:                     Computation effort in metric MADS

rp_K = rp/K;

D_ripples = (5.309e-03*log10(rp_K)^2 + 7.114e-02*log10(rp_K) - 0.4761)*log10(rc) - (2.66e-03*log10(rp_K)^2 + 0.5941*log10(rp_K) + 0.4278);

if K == 1,
    for i = 1 : K,
        alpha   = 2/(delta_F * prod(M(1:i)));
        beta(i) = M(i)/((prod(M(1:i))) * (1 - (((2-delta_F)/(2*OSR)) * prod(M(1:i)))));
        S       = alpha+sum(beta);
    end
else
    for i = 1 : K-1,
        alpha   = 2/(delta_F * prod(M(1:i)));
        beta(i) = M(i)/((prod(M(1:i))) * (1 - (((2-delta_F)/(2*OSR)) * prod(M(1:i)))));
        S       = alpha+sum(beta);
    end
end

RT = D_ripples * S * Fs;

if nargin > 7 & print_RT == 1,
    fprintf('The computational effort =%d\n', RT);
end