function [f, a, w, lengths] = MultiBandFilters(K, M, OSR, type)

%
% [f, a, w, lengths] = MultiBandFilters(K, M, OSR, type)
%
% This function design multi-band filters according to the input
% constarints. It supports two differnt types:
% 'NB' ->   In this type the filter design in the case of multi-stage will
%           be dependent on each stage, i.e., each stage will has its 
%           influence on the rest of the stages. Besides it has a narrow 
%           bandwidth  
%           Example:
%           K       = 3;
%           M       = [8 4 2];
%           OSR     = 64;
%           First stage bands  = [0 1/OSR 2/M(1)+-1/OSR ...]
%           Second stage bands = [0 1/OSR 2/M(1)*M(2)+-1/OSR ...]
%           .
%           .
%           and so on
%           As it can be observed the band here is 1/OSR
% 
% 'WB' ->   In this type the filter design in the case of multi-stage will
%           be independent on each stage, i.e., each stage will has 
%           no influence on the rest of the stages. 
%           Example:
%           K       = 3;
%           M       = [8 4 2];
%           OSR     = 64;
%           First stage bands  = [0 1/M(1) 2/M(1)+-1/M(1) ...]
%           Second stage bands = [0 1/M(2) 2/M(2)+-1/M(2) ...]
%           .
%           .
%           and so on
%           As it can be observed the band here is 1/M(i) where i is the
%           number of the stage, which means there is no influence of the
%           rest of the stages in each individual stage.
%
%       K:      Number od decimation stages
%       M:      Decimation stage factor for each individual stage
%       OSR:    Oversampling ratio
%       type :  'NB' Narrow Band
%               'WB' Wide Band
%
%       [f, a, w]: Exported parameters for remez|firpm Matlab function
%       lengths:   Defines the frequency band length for each filter stage,
%                  it is only used by the function named 'DecimationFilters'

if strcmp(type, 'NB'),
    
    for k = 1 : K-1,
        z = 0;
        y = 0;
        lengths(k) = prod(M(1:k))+2;
            for l = 1 : prod(M(1:k))+2,
                if l == 1,
                    f(k,l) = 0;
                    a(k,l)=1;
                else if l == 2
                    f(k,l) = 1/OSR;
                    a(k,l)=1;
                        else if mod(l,2)==1
                            f(k,l) = ((2+z)/prod(M(1:k)))-(1/OSR);
                            z = z + 2;
                                else if mod(l,2)==0
                                    f(k,l) = ((2+y)/prod(M(1:k)))+(1/OSR);
                                    y = y + 2;
                                        if l == prod(M(1:k))+2,
                                            f(k,l) = 1;
                                        end
                                            a(k,l)=0;
                                    end
                            end
                    end
                end
            end   
    end
    
elseif strcmp(type,'WB'),
    
    for k = 1 : K-1,
        z = 0;
        y = 0;
        lengths(k) = M(k)+2; %
            for l = 1 : M(k)+2; %
                if l == 1,
                    f(k,l) = 0;
                    a(k,l)=1;
                else if l == 2
                    f(k,l) = 1/(2*M(k));
                    %f(k,l) = 1/(M(k))-1/(2*OSR);
                    a(k,l)=1;
                        else if mod(l,2)==1
                            f(k,l) = (2+z)/M(k)-(1/(2*M(k)));
                            z = z + 2;
                                else if mod(l,2)==0
                                    f(k,l) = (2+y)/M(k)+(1/(2*M(k)));
                                    y = y + 2;
                                        if l == M(k)+2; %
                                            f(k,l) = 1;
                                        end
                                            a(k,l)=0;
                                    end
                            end
                    end
                end
            end   
    end
end
    
        
   c0 = 0.00000048003; % Added 00
   
   for m = 1 : K-1,
       x=0;
            for x = 1 : max(lengths)/2,
                wi(m,x) = c0/(2*sin(pi*f(m,2+x)/2))^3;
                x = x + 2;
            end
       wi(m,1)=0.00001; %% Changed 0.001
   end

    w = 1./wi;