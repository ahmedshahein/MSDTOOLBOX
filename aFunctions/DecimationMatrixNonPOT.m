function Dmatrix = DecimationMatrixNonPOT(D, K)

%
% deci_matrix = DecimationMatrixNonPOT(D, K)
%
% This function generate a matrix holding the possible decimation factors
% at each stage due to the overall decimation factor 'D' and number of
% decimation stages. The deciamtion factors are not constrained to be in a
% Power-Of-Two 'POT' format. It designed for maximum 4 decimation stages,
% i.e. K = [1 2 3 4].
%
%   D:  Decimation factor 'overall factor'
%   K:  Number of decimation stages
%   
%   deci_matrix:    Matrix of available decimation factor alternatives
%
% Example:
% D = 48, K = 3 -> deci_matrix = [[12 2 2]; [6 4 2]]

v = Factorize(D);
EvenOddFlag = EvenOdd(D);
deci_matrix = [];

if EvenOddFlag == 0,
    fprintf('WARNING: The decimation factor is an Odd number\n');
    fprintf('         The decimation factors might be not acceptable\n');
end

if isnan(v),
    fprintf('Non-acceptable decimation factor\n');
elseif length(v) <= 2,
    fprintf('WARNING: This decimation factor cannot hold more than 2 stages\n');
    Dmatrix = [D/2 2];
    if K == 1,
        deci_matrix = D;
    else 
        if CheckForInteger(D/2) == 1,
            deci_matrix = [D/2 2];
        else
            deci_matrix = v;
        end
    end    
else
    
    if K==3 & CheckForInteger(sqrt(D)) == 0 & (length(v) > (fix(log2(D)/2))),
        rows=fix(log2(D)/2);
    elseif K >= 3 & (length(v) == (fix(log2(D)/2)-1)),
        rows=fix(log2(D)/2)-2;
    else
        rows=fix(log2(D)/2)-1;
    end
    columns = K;
    matrix = zeros(rows, columns);
    matrix(:,1) = 2;

    if K == 1,
        deci_matrix = D;
    elseif K==2,
        deci_matrix = [2 D/2];
    elseif K == 3,    
       for j = 1 : rows,
             for i = 2 : K,
                if i == K,
                    matrix(j,i) = prod(v(i+j-1:length(v)));                
                    deci_matrix(j,:) = BubbleSort(matrix(j,:));
                else
                    matrix(j,i) = prod(v(i:i+j-1));
                    deci_matrix(j,:) = BubbleSort(matrix(j,:));
                end
            end
       end
    else     
       for j = 1 : rows,
             for i = 2 : K,
                if i == K,
                    matrix(j,i) = prod(v(i+j-1:length(v)));                
                    deci_matrix(j,:) = BubbleSort(matrix(j,:));
                elseif i == K-1,
                    matrix(j,i) = v(i+j-1);
                    deci_matrix(j,:) = BubbleSort(matrix(j,:));                
                else
                    matrix(j,i) = prod(v(i:i+j-1));
                    deci_matrix(j,:) = BubbleSort(matrix(j,:));
                end
            end
        end
    end 
    
    matrixd = fliplr(deci_matrix);
    
    if K > 2,
        if matrixd(rows,columns) == 1,
            Dmatrix = matrixd(1:rows-1,1:columns);
        else
            Dmatrix = matrixd;
        end
    else
        Dmatrix = matrixd;
    end

end


