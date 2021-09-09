function results = ex1()

% Example #1: We have a low Moes number and the fluid is isoviscous. This
% is the minimum amount of information that the solver needs to be able to
% produce a solution. At the end, we plot the resulting fluid pressure.
%
% To run, go to the project root directory and type: examples.ex1
%
% Copyright (C) 2021 theComputeKid

narginchk(0,0)

domain = setDomain();
moes = setMoes();

model = tribosolver(domain,moes);
results = model.solve();
results.plotP();

end

function domain = setDomain()

% Domain properties for the coarsest grid
nx = 32; xin = -5; xout = 1.25;
ny = 16; yin = -4; yout = 4;
mgl = 2; % Number of finer levels

domain = tribosolver.Domain(xin,xout,nx,yin,yout,ny,mgl);

end

function moes = setMoes()

% We set a low Moes M, with isoviscous fluid (L = 0)
M = 5; L = 0;

% Our initial solution guess of H0. A good guess can help solution
% stability.
H0 = -0.53;

moes = tribosolver.Moes(M,L,H0);

end
