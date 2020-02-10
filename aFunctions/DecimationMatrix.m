function m = DecimationMatrix(M,Stages)

%
%   m = DecimationMatrix(M,Stages)
%
% This function provide us with the differnet possibilties of recommended
% decimation stages according to number of decimation stages
%
%       M       = Input overall decimation factor
%       Stages  = Number of propossed decimation stages
%       matrix  = Output matrix, number of columns represent the number of
%                 stages of decimation, number of rows represent the number of different
%                 possibilities of decimation factors during the input decimation factor
%                 and number of stages
%
% Function practical specifications:
%   1. Decimation factors must be in descending factors
%   2. Input decimation factor must be power of 2
%   3. Last stage of decimation is always 2
%   4. Number of decimation stages are between 2-4
% The function is tested up to 2^15

rows = [];
columns = [];
p = 0;

% Single Stage decimation
if Stages == 1,
    m=M;
    
% Two stages decimation    
elseif Stages==2,
    m=[M/2,2];
    
% Three stages decimation 
else if Stages ==3 & rem(log2(M),2) == 0 | rem(log2(M),2) == 1,
    p=log2(M);
    rows=zeros;
    if mod(p,2)==1,
        rows=fix(p/2);
    else
        rows=p/2-1;
    end
    rows;
    columns=Stages;
    matrix=zeros(rows,columns);
    
    for j=1:rows,
        matrix(j,3)=2;
    end
    
    n=1;
    for j=1:rows,
        matrix(j,2)=2^n;
        n = n+1;
    end   
    
    if mod(p,2)==0,
        k=rows-1;
        for j=1:rows,
            matrix(j,1)=sqrt(M)*2^k;
            k = k-1;
        end
       
    else
        k=rows-1;
        for j=1:rows,
            matrix(j,1)=sqrt(M/2)*2^k;
            k = k-1;
        end
           
    end
    m = matrix;

% Four stages decimation
else if Stages == 4 & rem(log2(M),2) == 0 | rem(log2(M),2) == 1,
    columns = Stages;
    rows = zeros;
    p=log2(M);
    
    if p == 4,
        rows = 1;
    elseif OrCondition(5,9,p)==1,
        rows = p-4;
    elseif OrCondition(10,11,p)==1,
        rows = p-3;
    elseif p == 12,
        rows = p-2;
    elseif OrCondition(13,15,p)==1,
        if mod(ceil(p/2),2)==0,
            rows = ceil(p/2)*2;
        else
            rows = fix(p/2)*2;
        end
    end
end

matrix=zeros(rows,columns);
    if p==4,
        m=[2,2,2,2];
    else
    
    for j=1:rows,
        matrix(j,4)=2;
    end
    
    h=1;
    k=fix(p/2)-1;
        if mod(p,2)==0,
            for j = 1:round(p/2)-1,
                matrix(j,3)=2;
                matrix(j,2)=2^h;
                h=h+1;
                matrix(j,1)=(sqrt(M)/2)*2^(k-1);
                k = k-1;
            end
            
            h=2;
            k=fix(p/2)-1;
            for j = round(p/2):(round(p/2)+(k-3)),
                matrix(j,3)=4;
                matrix(j,2)=2^h;
                h=h+1;
                matrix(j,1)=(sqrt(M)/8)*2^(k-1);
                k = k-1;
            end
           
            h=3;
            k=fix(p/2)-1;
            for j = (round(p/2)+(k-2)):(round(p/2)+(k-2)+(k-4)),
                matrix(j,3)=8;
                matrix(j,2)=2^h;
                h=h+1;
                matrix(j,1)=(sqrt(M)/32)*2^(k-1);
                k = k-1;
            end
            
            if p>13,
                h=4;
                matrix(rows,3)=16;
                matrix(rows,2)=2^h;
                h=h+1;
                matrix(rows,1)=(sqrt(M)/32)*2^(k);
            end
          
        else
               for j = 1:fix(p/2)-1,
                matrix(j,3)=2;
                matrix(j,2)=2^h;
                h=h+1;
                matrix(j,1)=(sqrt(M/2)/2)*2^(k);
                k = k-1;
            end
    
     
            h=2;
            k=fix(p/2)-1;
            for j = fix(p/2):(fix(p/2)+(k-2)),
                matrix(j,3)=4;
                matrix(j,2)=2^h;
                h=h+1;
                matrix(j,1)=(sqrt(M/2)/8)*2^(k);
                k = k-1;
            end
          
            
            h=3;
            k=fix(p/2)-1;
            for j = (fix(p/2)+(k-1)):(fix(p/2)+(k-1)+(k-4)),
                matrix(j,3)=8;
                matrix(j,2)=2^h;
                h=h+1;
                matrix(j,1)=(sqrt(M/2)/16)*2^(k-1);
                k = k-1;
            end
           
     
            if p>12,
            h=4;
            k=fix(p/2)-1;
            for j = (fix(p/2)+(k-1)+(k-3)):rows,
                matrix(j,3)=16;
                matrix(j,2)=2^h;
                h=h+1;
                matrix(j,1)=(sqrt(M/2)/64)*2^(k-1);
                k = k-1;
            end
        end
        
        end
        m = matrix;
      end
  end
end