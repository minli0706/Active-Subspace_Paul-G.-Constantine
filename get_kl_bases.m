function [U,sv] = get_kl_bases(corr_length,m,pde_data,filename)

% corr_length: correlation length
% m: truncation of KL
% pde_data: mesh for PDE

filename = sprintf('kl/%s',filename);

% Random field model for diffusion coefficients.
if exist(filename,'file')
    % If the PDEs have been solved, use 'em.
    load(filename);
else
    corr.name='exp';
    corr.c0=[corr_length corr_length];
    corr.sigma=1;
    [~,KL]=randomfield(corr,pde_data.p','trunc',m);
    U=bsxfun(@times,KL.bases,KL.sv');
    sv = KL.sv;
    save(filename,'U','sv');
end