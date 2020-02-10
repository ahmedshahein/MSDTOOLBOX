function HB_stage = HalfBandFilter(Rpb, Rsb, Type, deltaF_2, N)

%
% HB_stage = HalfBandFilter(Rpb, Rsb, Type, deltaF_2, N)
%
% This function is used to design the half-band filter in 2 different
% algorithms.
% 
% Rpb   : Pass-band ripples
% Rsb   : Stop-band attenuation
% Type  : Half-band type I or II
% N     : Filter order
%
% P.S.
% For Type-I, it is prefered to use the ripples in dB format not linear
% format.

% Assure that the ripples are in linear format.
if Rsb > 1,
    [rp rs] = Ripples(Rpb, Rsb, 'dB2L');
else
    rp = Rpb;
    rs = Rsb;
end

% Designing Half-band stage
if Type == 'Type-I ',
    Fhb     = [0.5-deltaF_2 0.5+deltaF_2]; 
    ahb     = [1 0];
    devhb   = [rp rs]; 
    [nhb,fohb,aohb,whb] = remezord(Fhb, ahb, devhb);
    bhb = remez(nhb,fohb,aohb,whb);
    HB_stage = ZeroOutCoeff(bhb);
else
    if mod(N,2) == 0,
        HB_stage = 'WARNING: The filter order has to be ODD!!!'
        return
    end
    if isempty(N),
        HB_stage = 'WARNING: Enter the filter order!!!'
        return
    else
        P = (N+1)/2; 
        Filter_Coeff_Length = [1:N];
        Filter_Coeff_Length = Filter_Coeff_Length - P;
        HB_Sinc_Filter = sin(pi*Filter_Coeff_Length/2)./(pi*Filter_Coeff_Length);
        HB_Sinc_Filter(P)=.5;
        thewin=window(@hamming, N);
        HB_Window = HB_Sinc_Filter .* thewin';
        HB_Filter = HB_Window/sum(HB_Window);
        HB_stage = zerooutcoeff(HB_Filter);
    end
end