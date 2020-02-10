function  adders = NumberOfAddersInMuliplier(coeff, n)

%
% adders = NumberOfAddersInMuliplier(coeff, n)
%
% Calculate number of adders in a whole multiplier according to the number
% of binary multipliers. This number is from the hardware point of view,
% because at each multiplier only ones will be implemented as hardwired
% multiplier. We use it to estimate the number of adders for FIR filter
%
%   coeff:    Normalized filter coefficients
%   n:        Quantization bitwidth for the declared coefficients
%
%   adders:   Number of adders in constant multiplier

bin_array = dec2bin_array_2s(coeff, n);
 
num_binary_multipliers = sum(bin_array);

num_adders = num_binary_multipliers - 1;

adders = sum(num_adders);

%%%%%%%%%%%%%%%
% SUBFUNCTION %
%%%%%%%%%%%%%%%
function bin_array = dec2bin_array_2s(dec_array, n)

% This function convert decimal vector to its equivelent binary array, and
% it support negative numbers too represented in 2's complement format

% dec_array :   input decimal vector, each element in the vector represent
%               single decimal number
% n         :   binary bit width
% bin_array :   output binary array equivelent to the input decimal vecctorfunction bin_array = dec2bin_array_2s(dec_array, n);

bin_array = zeros(length(dec_array), n);

for i = 1 : length(dec_array),
    bin_array(i,:) = decimal2binary_2s(dec_array(i), n);
end

%%%%%%%%%%%%%%%
% SUBFUNCTION %
%%%%%%%%%%%%%%%
function binary = decimal2binary_2s(decimal, n)

% This function convert decimal number to its 2's complement binary equivelent

tmp = decimal;
count = 0;

for i = 1 : n,
    if i == n,
        one(i) = 1;
    else 
        one(i) = 0;
    end
end

while tmp > 1,
    tmp = tmp/2;
    count = count + 1;
end
min_length = count;

if n < min_length,
    fprintf('The number of bits is not sufficent \n');
else
    if decimal >= 0,
        dividend = decimal;
        for i = 1 : n,
            remainder(i) = rem(dividend, 2);
            dividend = floor(dividend/2);
        end
            binary = fliplr(remainder);
    else
        dividend = abs(decimal);
        for i = 1 : n,
            remainder(i) = rem(dividend, 2);
            dividend = floor(dividend/2);
        end
            norm_binary = fliplr(remainder);  
            complement = not(norm_binary);
            binary_tmp = complement + one;
            
            for i = 1 : length(binary_tmp),
                if binary_tmp(i) > 1,
                    flag = 1;
                else
                    flag = 0;
                end
            end
            
            if flag == 1,
                binarytmp = binary2decimal(complement) + 1;
                binary = decimal2binary(binarytmp, n);
            else
                binary = complement + one;
            end
    end
end

%%%%%%%%%%%%%%%
% SUBFUNCTION %
%%%%%%%%%%%%%%%
function binary = decimal2binary(decimal, n)

% This function convert decimal number to its binary equivelent

dividend = decimal;
tmp = decimal;
count = 0;

while tmp > 1,
    tmp = tmp/2;
    count = count + 1;
end
min_length = count;

if n < min_length,
    fprintf('The number of bits is not sufficent \n');
else
    for i = 1 : n,
        remainder(i) = rem(dividend, 2);
        dividend = floor(dividend/2);
    end
        binary = fliplr(remainder);
end

%%%%%%%%%%%%%%%
% SUBFUNCTION %
%%%%%%%%%%%%%%%
function decimal = binary2decimal(binary);

% This function binary number to its decimal equivelent
% binary    : input binary number with any bit width
% decimal   : output decimal equivelent to the binary input

binary = fliplr(binary);

for i = 1 : length(binary),
    if binary(i) == 1,
        tmp(i) = 2^(i-1);
    else
        tmp(i) = 0;
    end
    
end
decimal = sum(tmp);