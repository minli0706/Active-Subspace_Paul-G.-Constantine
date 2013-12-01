function dnegloglik = d_neg_log_likelihood(X,y,lambda,gamma)

dtol = 1e-6;
dnegloglik = (0.5/dtol)*(neg_log_likelihood(X,y,lambda,gamma+dtol)...
    -neg_log_likelihood(X,y,lambda,gamma-dtol));