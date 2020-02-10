function [FilteredData Hcic] = CICFilter(M, N, R, Fs, Data)

%
% [FilteredData Hcic] = CICFilter(M, N, R, Fs, Data)
%
% This function formulate the CIC transfer function, which can be used to
% filter the filter input data and plot the filter response with both input
% and output spectrums. It is integrated with the 
% 'CascadedCICCompensationFilter' function.
%
%   M:              Differntial delay
%   N:              Filter order
%   R:              Decimation factor
%   Fs:             Sampling frequency
%   Data:           Filter input data
%
%   FilteredData:   Filtered and downsampled output data
%   Hcic:           CIC filter transfer function
%
% Ref: GUI Based Decimation Filter Design Toolbox for Multi Standard Wireless Transceveris

Num     = ones(1,R*M);
Den     = 1;

H       = filt(Num,Den,1/Fs);

Hcic    = (H^N)/((R*M)^N);

Num_CIC = cell2mat(Hcic.num);
Den_CIC = cell2mat(Hcic.den);

FilteredData = downsample(filter(Num_CIC,Den_CIC,Data),R);