function results = ex3()

% Example #3: Non-Isoviscous fluids with graphical verbosity (caution: This
% plotting slows down the solution process)
%
% To run, go to the project root directory and type: examples.ex3
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

nx = 32; xin = -3; xout = 1.5;
ny = 32; yin = -2.5; yout = 2.5;
mgl = 2;

domain = tribosolver.Domain(xin,xout,nx,yin,yout,ny,mgl);

end

function moes = setMoes()

M = 15; L = 5;
H0 = -0.53;

moes = tribosolver.Moes(M,L,H0);

end

function exec = setExecution()

% We solve using double precision using the CPU
prec = "double"; dev = "cpu";

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
gsSORFactor = 3e-1;
h0UpdateFactor = 5e-3;

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