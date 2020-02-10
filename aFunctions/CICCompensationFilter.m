function Hcomp = CICCompensationFilter(N, R, Fs)

% 
% Hcomp = CICCompensationFilter(N, R, Fs)
%
% This function formulate the CIC compensation filter transfer function.
% It is integrated with the 'CascadedCICCompensationFilter' function.
%
%   N:              Filter order
%   R:              Decimation factor
%   Fs:             Sampling frequency
%
%   Hcomp:          CIC compensation filter transfer function
%
% Ref: Design of CIC Compensator Filter in a Digital IF Receiver

b = [2 1 0 0 0 -1 -2];
b = b(N);

A = -1 * (2^(b+2) + 2);
B = -2^(-1*(b + 2));

Num = zeros(1,2*R+1);
Num(1) = 1;
Num(R+1) = A;
Num(end) = 1;

Hcomp = filt(B*Num,1,1/Fs);

figure
PlotFreqResponse(cell2mat(Hcomp.num), Fs)