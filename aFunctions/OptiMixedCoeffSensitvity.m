function [test_filter_coeff_sensitvity count_indices]= OptiMixedCoeffSensitvity(coeff, Q, Sn, Sn_value, N_NPT, Fs)

filter_coeff_indices            = [];
count_indices                   = 0;

% Defining the location of the coefficients to be rounded to NPT factors
% for each sensitvity factor, and the number of coefficients to be changed.
if Sn_value == 0,

    Sn_elements = ElementsOfVector(Sn);

    for i = 1 : length(Sn_elements), 
        Sn_max = Sn_elements(i);
        count = 0;
        for j = 1 : length(coeff),
            if Sn(j)>0 & Sn(j)<=Sn_max,
                filter_coeff_indices(i,j) = j;
                count_indices(i) = count + 1;
            else
                filter_coeff_indices(i,j) = 0;
            end
        end
    end

else

    Sn_max = Sn_value;
    
    for k = 1 : length(coeff),
        if Sn(k)>0 & Sn(k)<=Sn_max,
            filter_coeff_indices(k) = k;
            count_indices = count_indices + 1;
        else
            filter_coeff_indices(k) = 0;
        end
    end
end

% Rounding the low sensitive coefficients to its NPT factor, whether for N =
% 1, 2, or 3.
if Sn_value ~= 0,
    
    for l = 1 : length(coeff),
        if N_NPT == 1,
            if filter_coeff_indices(l) ~= 0,
                test_filter_coeff_sensitvity(l)  = SPT(NormalizedCoeff(coeff, Q, coeff(l)));
            else
                test_filter_coeff_sensitvity(l)  = NormalizedCoeff(coeff, Q, coeff(l));
            end
        elseif N_NPT == 2,
            if filter_coeff_indices(l) ~= 0,
                test_filter_coeff_sensitvity(l)  = NPT(NormalizedCoeff(coeff, Q, coeff(l)),N_NPT,Q);
            else
                test_filter_coeff_sensitvity(l)  = NormalizedCoeff(coeff, Q, coeff(l));
            end
        elseif N_NPT == 3,
            if filter_coeff_indices(l) ~= 0,
                test_filter_coeff_sensitvity(l)  = NPT(NormalizedCoeff(coeff, Q, coeff(l)),N_NPT,Q);
            else
                test_filter_coeff_sensitvity(l)  = NormalizedCoeff(coeff, Q, coeff(l));
            end
        end        
    end
    
else

    [rows columns] = size(filter_coeff_indices);
    
    for m = 1 : rows,
        for n = 1 : columns,
            if N_NPT == 1,
                if filter_coeff_indices(m,n) ~= 0,
                    test_filter_coeff_sensitvity(m,n)  = SPT(NormalizedCoeff(coeff, Q, coeff(n)));
                else
                    test_filter_coeff_sensitvity(m,n)  = NormalizedCoeff(coeff, Q, coeff(n));
                end
            elseif N_NPT == 2,
                if filter_coeff_indices(m,n) ~= 0,
                    test_filter_coeff_sensitvity(m,n)  = NPT(NormalizedCoeff(coeff, Q, coeff(n)),N_NPT,Q);
                else
                    test_filter_coeff_sensitvity(m,n)  = NormalizedCoeff(coeff, Q, coeff(n));
                end
            elseif N_NPT == 3,
                if filter_coeff_indices(m,n) ~= 0,
                    test_filter_coeff_sensitvity(m,n)  = NPT(NormalizedCoeff(coeff, Q, coeff(n)),N_NPT,Q);
                else
                    test_filter_coeff_sensitvity(m,n)  = NormalizedCoeff(coeff, Q, coeff(n));
                end
            end             
        end
    end
    
end

% FIG = figure('Name', 'Sensitvity of Filter Coefficients', 'NumberTitle' , 'off');
% plot(abs(Sn));
% grid on
% xlabel('Filter Coefficeints');
% ylabel('Coefficients Sensitvity');
% title('Sensitvity of Filter Coefficients');
% 
% FIG = figure('Name', 'Filters Frequency Responses', 'NumberTitle' , 'off');
% filters = test_filter_coeff_sensitvity;
% plot_freq_response(filters, Fs);
% title('Filter Response');
% legend('Normalized Coeff.','NPT Coeff.','Mixed Coeff.');

% End