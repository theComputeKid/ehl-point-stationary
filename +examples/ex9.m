function results = ex9()

% Example #9: Uses the gpu parallel line solver with multiple grids and
% visual feedback.
%
% To run, go to the project root directory and type: examples.ex9
%
% Copyright (C) 2021 theComputeKid

narginchk(0,0)

domain = setDomain();
moes = setMoes();
exec = setExecution();
relax = setRelaxation();

% The command to initialize the solver can take arguments in any order, so
% compared to Example 2, we switch around the places of the inputs to prove
% our point.
model = tribosolver(relax,domain,exec,moes);
results = model.solve();
results.plotP();

end

function domain = setDomain()

nx = 64; xin = -2.5; xout = 1.5;
ny = 64; yin = -2.5; yout = 2.5;
mgl = 5;

domain = tribosolver.Domain(xin,xout,nx,yin,yout,ny,mgl);

end

function moes = setMoes()

M = 15; L = 4;
H0 = -0.63;

moes = tribosolver.Moes(M,L,H0);

end

function exec = setExecution()

% We solve using single precision using the GPU sequential solver.
prec = "single"; dev = "gpu_par";

% A verbosity level of 2 indicates the display of both text (verbosity > 0)
% and graphical (verbosity > 1) plots during the solution scheme. Note that
% the display of graphical plots is computationally expensive.
verbosity = 2;

exec = tribosolver.Execution(prec,dev,verbosity);
end

function relax = setRelaxation()

% The solution method uses both Gauss-Seidel and Jacobian iterations. Here
% we specify the successive-over-relaxation paramters for each type of
% update.
jacobianSORFactor = 6e-1;
gsSORFactor = 8e-1;
h0UpdateFactor = 1e-2;
epsSwitch = 0.3;

% We set a maximum number of FMG Cycles (in-case our solution does not
% converge)
numCycles = 500;

% A FMG V-Cycle would be gamma=1, while an FMG W-Cycle would be gamma=2
gamma = 2;

% Number of MultiGrid Iterations
itPre = 5;
itPost = 5;
itMain = 10;

relax = tribosolver.Relaxation( ...
    jacobianSORFactor,gsSORFactor,h0UpdateFactor, ...
    numCycles,gamma, ...
    epsSwitch, ...
    itPre, itMain, itPost ...
    );

end