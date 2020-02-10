function o = VectorContentIterations(x)

%
% o = VectorContentIterations(x)
%
% This function exports statistical information about input vector x.
%   Example:
%   x = [5 8 3 4 9 13 4];
%   o = vector_content_iterations(x);
%   o = 
%       3     1
%       4     2
%       5     1
%       8     1
%       9     1
%      13     1
    
outv = zeros(max(x)-min(x)+1,2);
for i=min(x):max(x)
    outv(i,1) = i;
    outv(i,2) = length(find(x==i));
end

exact=0;
for i =1:length(outv),
    if outv(i,2)~=0,
        exact=exact+1;
    end
end

o=zeros(exact,2);
k=1;
for j = 1:length(outv),
    if outv(j,2)~=0,
        o(k,:)=outv(j,:);
        k=k+1;
    end
end
    
    