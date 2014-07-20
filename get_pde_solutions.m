function [f,G] = get_pde_solutions(X,U,pde_data,filename)

if exist(sprintf('pde/%s',filename),'file') && ~isempty(filename)
    % If the PDEs have been solved, use 'em.
    load(sprintf('pde/%s',filename));
else
    
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

    f=zeros(M,1);

    if nargout>1, G=zeros(M,m); else, G=[]; end
    for i=1:M
        tic

        % forward
        c=pdeintrp(p,t,exp(U*X(i,:)'));
        [K,F]=assempde(b,p,e,t,c,0,1);

        % quantity of interest
        u = K\F;
        f(i)=psi'*u;

        if nargout>1
            % Get gradients of objective function

            % adjoint
            phi=K\psi;

            % sensitivities
            C=bsxfun(@times,exp(U*X(i,:)'),U);
            for j=1:m
                [K,~]=assempde(b,p,e,t,...
                    pdeintrp(p,t,C(:,j)),0,1);
                G(i,j)=-phi'*(K*u);
            end
        end

        fprintf('PDE solves: %i of %i, time: %4.2f\n',i,M,toc);

    end
    
    if ~isempty(filename), save(sprintf('pde/%s',filename),'G','f','psi','X'); end

end

    