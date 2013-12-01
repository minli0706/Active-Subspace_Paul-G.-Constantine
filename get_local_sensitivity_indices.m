function ind = get_local_sensitivity_indices(X,U,pde_data)

p = pde_data.p;
e = pde_data.e;
t = pde_data.t;
b = pde_data.b;
[M,m] = size(X);

% defining the QOI as linear functional of solution.

psi=zeros(size(p,2),1);
psi(p(1,:)==1)=1; 
[~,A,~]=assema(p,t,0,1,0);
psi = A*psi;

z = zeros(m,1);

% forward
c=pdeintrp(p,t,exp(U*z));
[K,F]=assempde(b,p,e,t,c,0,1);

% quantity of interest
u = K\F;

% adjoint
phi=K\psi;

% sensitivities
g = zeros(m,1); 
C=bsxfun(@times,exp(U*z),U);
for j=1:m
    [K,~]=assempde(b,p,e,t,...
        pdeintrp(p,t,C(:,j)),0,1);
    g(j)=-phi'*(K*u);
end

[~,ii] = sort(abs(g),'descend');
ind = ii(1:2);
