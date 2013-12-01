function [gstar,wstar,parms] = gpml_regression(X,y,Xstar,gamma_bnds,lambda)

% N is number of samples
% n is dimension of subspace
[N,n] = size(X);

if nargin==3
    % wrapper for the gmpl code

    cov = {@covSEiso}; sf = 1; ell = 2;
    hyp0.cov  = log([ell;sf]);
    mean = {@meanSum,{{@meanPoly,2},@meanConst}};
    if n==1
        mean = {@meanSum,{{@meanPoly,2},@meanConst}};
        hyp0.mean = [0; 0; 0];
    elseif n==2
        mean = {@meanSum,{{@meanPoly,2},@meanConst}};
        hyp0.mean = [0; 0; 0; 0; 0];
    else
        mean = {@meanSum,{{@meanPoly,1},@meanConst}};
        hyp0.mean = zeros(n+1,1);
    end

    Ncg = 500;
    lik = 'likGauss';
    hyp0.lik  = log(0.2);
    inf = 'infExact';

    hyp = minimize(hyp0,'gp', -Ncg, inf, mean, cov, lik, X, y);
    [gstar, wstar] = gp(hyp, inf, mean, cov, lik, X, y, Xstar);
    negloglik = gp(hyp, inf, mean, cov, lik, X, y);
    fprintf('Negative log likelihood: %6.4f\n',negloglik);
    parms=hyp;
    
else
    % Custom code for heuristics
    
    fun = @(g) neg_log_likelihood(X,y,lambda,g);
    gamma = fminbnd(fun,gamma_bnds(1),gamma_bnds(2),...
        optimset('Display','iter','MaxFunEvals',500,'TolX',1e-5));
    
    % bisection with finite differences for minimization
    %a=gamma_bnds(1); b=gamma_bnds(2);
    %fun=@(g) d_neg_log_likelihood(X,y,lambda,g);
    %btol = 1e-5; funtol=1e-5; count=1; maxcount=500;
    %fprintf('Bisection bounds: a=%6.4f, b=%6.4f\n',a,b);
    %while(abs(b-a) >= btol || (abs(fun(a)) >= funtol && abs(fun(b)) >=funtol))
    %    c=0.5*(a+b);
    %    if fun(c) == 0;
    %        break;
    %    elseif fun(a)*fun(c) > 0
    %        b=c;
    %    else
    %        a=c;
    %    end
    %    fprintf('%0.4d, neg-log-lik: %6.4f, c: %6.4f, fun(a): %6.4f, fun(b): %6.4f\n',...
    %        count,neg_log_likelihood(X,y,lambda,c),c,fun(a),fun(b));
    %    count=count+1;
    %    if count>maxcount, break; end
    %end
    %gamma = c;
    
    %ngamma = 20;
    %gammas = linspace(gamma_bnds(1),gamma_bnds(2),ngamma);
    %neglogliks = zeros(ngamma,1);
    %for i=1:ngamma
    %    sigma = gammas(i)*sum(lambda);
    %    sigma_n = gammas(i)*sum(lambda(n+1:end));
    %    ell = sigma./lambda(1:n);
    %    neglogliks(i) = neg_log_likelihood(X,y,sigma_n,ell,sigma);
    %    fprintf('Negative log likelihood(%d): %6.4f\n',i,neglogliks(i));
    %end
    %[~,ind] = min(neglogliks);
    %gamma = gammas(ind);
    
    sigma = gamma*sum(lambda);
    sigma_n = gamma*sum(lambda(n+1:end));
    ell = sigma./lambda(1:n);
    parms.gamma=gamma; parms.sigma=sigma;
    parms.sigma_n=sigma_n; parms.ell=ell;

    % covariance matrix of observations
    Ky = exp2_cov(X,X,sigma,ell) + sigma_n*eye(N);

    % covariance of test sites
    kstar = exp2_cov(X,Xstar,sigma,ell);

    % Following algorithm 2.1 in Rasmussen
    L = chol(Ky,'lower'); % cholesky
    alpha = L'\(L\y);

    % prediction and variance without polys
    fstar = kstar'*alpha;
    vstar = diag(exp2_cov(Xstar,Xstar,sigma,ell)) - sum(((L\kstar)').^2,2);

    % coefficients of polynomial bases
    H = poly2_bases(X);
    A = H'*(Ky\H);
    beta = A\(H'*(Ky\y));

    % prediction and variance with polys
    Hstar = poly2_bases(Xstar);
    R = Hstar' - H'*(Ky\kstar);
    gstar = fstar + R'*beta;
    wstar = vstar + diag(R'*(A\R));
end