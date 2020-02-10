function [DataOut Hcascaded] = CascadedCICCompensationFilter(M, R, N, Fs, Data)

%
% [DataOut Hcascaded] = CascadedCICCompensationFilter(M, R, N, Fs, Data)
%
% This function design a CIC filter cascaded with compensation filter.
% The design parameters is given by:
%
%   M:              Differntial delay
%   N:              Filter order
%   R:              Decimation factor
%   Fs:             Sampling frequency
%   Data:           Filter input data
%
%   DataOut:        Filtered input data (Data) using the the 
%                   design CIC filter
%   Hcascaded:      Transfer function of cascaded CIC (CIC+Compensation)

[FilteredData Hcic] = CICFilter(M, N, R, Fs, Data);

Hcomp = CICCompensationFilter(N, R, Fs);

Hcascaded = Hcic * Hcomp;

DataOut = downsample(filter(cell2mat(Hcascaded.num), cell2mat(Hcascaded.den), Data), R);

figure
PlotFreqResponse(cell2mat(Hcascaded.num), Fs,'k')