function results = ex2()

% Example #2: Same as Example 1 except:
%
% 1. We choose a verbosity level that shows textual results during the
% solution process.
%
% 2. We provide our own relaxation factors
%
% 3. We plot the film thickness from the results
%
% To run, go to the project root directory and type: examples.ex2
%
% Copyright (C) 2021 theComputeKid

narginchk(0,0)

domain = setDomain();
moes = setMoes();
exec = setExecution();
relax = setRelaxation();

model = tribosolver(domain,moes,exec,relax);
results = model.solve();
results.plotH();

end

function domain = setDomain()

nx = 32; xin = -4.5; xout = 1.5;
ny = 16; yin = -3; yout = 3;
mgl = 3;

domain = tribosolver.Domain(xin,xout,nx,yin,yout,ny,mgl);

end

function moes = setMoes()

% We set a low Moes M, with isoviscous fluid (L = 0)
M = 10; L = 0;

% Our initial solution guess of H0. A good guess can help solution
% stability.
H0 = -0.53;

moes = tribosolver.Moes(M,L,H0);

end

function exec = setExecution()

% We solve using double precision using the CPU
prec = "double"; dev = "cpu_seq";

% A verbosity level of 2 indicates the display of both text (verbosity > 0)
% and graphical (verbosity > 1) plots during the solution scheme. Note that
% the display of graphical plots is computationally expensive.
verbosity = 1;

exec = tribosolver.Execution(prec,dev,verbosity);
end

function relax = setRelaxation()

% The solution method uses both Gauss-Seidel and Jacobian iterations. Here
% we specify the successive-over-relaxation paramters for each type of
% update.
jacobianSORFactor = 6e-1;
gsSORFactor = 8e-1;
h0UpdateFactor = 5e-2;

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