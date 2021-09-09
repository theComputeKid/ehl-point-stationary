function results = ex4()

% Example #4: Non-Isoviscous fluids in single precision. Setting many other
% custom properties that the solver exposes. Caution: Slightly long-running
% example (may take > 5 mins).
%
% To run, go to the project root directory and type: examples.ex4
%
% Copyright (C) 2021 theComputeKid

narginchk(0,0)

domain = setDomain();
moes = setMoes();
exec = setExecution();
relax = setRelaxation();
tol = setTolerance();

model = tribosolver(exec,relax,tol,domain,moes);
results = model.solve();
results.plotP();

end

function domain = setDomain()

% If a dense grid is required to solve a particular problem, prefer
% increasing the density in X (nx) rather than Y (ny). This is because the
% solver solves one line of Y at a time in a vectorized manner, which is
% then looped over all the lines. The solver is more efficient at solving a
% single line (i.e. for all points 1:nx), than looping over all lines
% (1:ny).
nx = 32; xin = -3; xout = 1.25;
ny = 16; yin = -2.25; yout = 2.25;
mgl = 4;

domain = tribosolver.Domain(xin,xout,nx,yin,yout,ny,mgl);

end

function moes = setMoes()

M = 20; L = 5;
H0 = -0.8;

moes = tribosolver.Moes(M,L,H0);

end

function exec = setExecution()

% We solve using single precision using the CPU
prec = "single"; dev = "cpu";

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
jacobianSORFactor = 2e-1;
gsSORFactor = 2e-1;
h0UpdateFactor = 5e-4;

% We set a maximum number of FMG Cycles (in-case our solution does not
% converge)
numCycles = 50;

% A FMG V-Cycle would be gamma=1, while an FMG W-Cycle would be gamma=2
gamma = 2;

relax = tribosolver.Relaxation( ...
    jacobianSORFactor,gsSORFactor,h0UpdateFactor, ...
    numCycles,gamma ...
    );

end

function tol = setTolerance()

% Set a custom tolerance for each solution variable
p = 1e-2;
fb = 1e-2;
h = 1e-2;

tol = tribosolver.Tolerance(p,h,fb);

end