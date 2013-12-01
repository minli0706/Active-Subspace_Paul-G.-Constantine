function negloglik = neg_log_likelihood(X,y,lambda,gamma)

[M,n] = size(X);
sigma = gamma*sum(lambda);
sigma_n = gamma*sum(lambda(n+1:end));
ell = sigma./lambda(1:n);

% covariance matrix of observations
Ky = exp2_cov(X,X,sigma,ell) + sigma_n*eye(M);
L = chol(Ky,'lower');

% polynomial bases
H = poly2_bases(X);
A = H'*(Ky\H);
B = chol(A,'lower');

% negative log likelihood
z = Ky\y;
negloglik = 0.5*(y'*(Ky\y)) ...
    - 0.5*((H'*z)'*(A\(H'*z))) ...
    + 0.5*sum(log(diag(L))) ...
    + 0.5*sum(log(diag(B))) ...
    + 0.5*(M-size(H,2))*log(2*pi);
    