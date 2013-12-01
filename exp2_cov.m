function C = exp2_cov(X,Y,sigma,ell)

m=size(X,1); n=size(Y,1);

if isscalar(ell), ell=ell*ones(size(X,2),1); end
c=-0.5./ell; c=c(:);
C=zeros(m,n);
for i=1:n
    point=Y(i,:);
    XX=bsxfun(@minus,X,point);
    C(:,i)=sigma.*exp((XX.^2)*c);
end