function [rpb rsb] = Ripples(rp, rc, convert)

%
% [rpb rsb] = Ripples(rp, rc, convert)
%
% This function to convert the passband ripples and stopbanb attenuation 
% from or to Linear format to or from dB format.
%
%   rp:         Passband ripples        'absoulte value without signs'
%   rc:         Stopband attenuation    'absoulte value without signs'
%   convert:    Flag to convert inputs form 
%               'L2dB'  Linear -> dB    
%               'dB2L'  dB -> Linear
%

if convert == 'L2dB',
    rp_dB = 20*log10((1+rp)/(1-rp));
    rc_dB = -20*log10(rc);
    rpb = rp_dB;
    rsb = rc_dB;
elseif convert == 'dB2L',
    rp_L=(10^(rp/20)-1)/(10^(rp/20)+1);
    rc_L=10^(-rc/20);
    rpb = rp_L;
    rsb = rc_L;
end