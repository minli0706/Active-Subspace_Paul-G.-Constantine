function data = get_pde_data()
% Returns PDE geometry and mesh specs in a struct

% load geometry and boundary conditions
bg=load('bg.mat'); b=bg.b; g=bg.g; clear bg;
data.g = g; data.b = b;

% Square physical domain.
[p,e,t]=initmesh(g,'hmax',0.01);
data.p = p;
data.e = e;
data.t = t;