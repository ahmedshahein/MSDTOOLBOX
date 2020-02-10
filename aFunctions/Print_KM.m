function Print_KM(K, M)

if K == 2,
    fprintf('\nOptimal number of decimation stages K = %d', K);
    fprintf('\nOptimal decimation factor for each deciamtion stage = [%d %d]\n', M);
elseif K == 3,
    fprintf('\nOptimal number of decimation stages K = %d', K);
    fprintf('\nOptimal decimation factor for each deciamtion stage = [%d %d %d]\n', M);    
elseif K == 4,
    fprintf('\nOptimal number of decimation stages K = %d', K);
    fprintf('\nOptimal decimation factor for each deciamtion stage = [%d %d %d %d]\n', M);    
end
