function H = poly2_bases(X)

[M,d] = size(X);

% build polynomial constraints
I=index_set('full',2,d); 
I = flipud(I);
nn=size(I,2);
H=zeros(M,nn);
for i=1:nn
    H(:,i)=prod(bsxfun(@power,X,I(:,i)'),2);
end