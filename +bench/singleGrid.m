function singleGrid()

% Benchmarks
%
% To run, go to the project root directory and type: bench.singleGrid
%
% Copyright (C) 2021 theComputeKid

narginchk(0,0)

domain = setDomain();
moes = setMoes();
relax = setRelaxation();
tol = setTolerance();
disp("------------")
disp("SINGLE GRID BENCHMARKS:")
disp("Total Cycles: " + relax.numCycles)
disp("Coarsest Grid Density: " + domain.nx + " x " + domain.ny)

for i = ["cpu_par", "cpu_seq", "gpu_par"]
    disp("------------")
    disp("CASE: " + i)
    disp("------------")
    exec = setExecution(i);
    model = tribosolver(domain,moes,exec,relax,tol);
    model.solve();
end

end

function tol = setTolerance()
p = 5e-3; h = 1e-2; fb = 1e-2;
tol = tribosolver.Tolerance(p,h,fb);
end

function domain = setDomain()

nx = 1024; xin = -4.5; xout = 1.5;
ny = 1024; yin = -3; yout = 3;
mgl = 1;

domain = tribosolver.Domain(xin,xout,nx,yin,yout,ny,mgl);

end

function moes = setMoes()
M = 10; L = 2;
H0 = -0.53;
moes = tribosolver.Moes(M,L,H0);
end

function exec = setExecution(dev)
prec = "single";
verbosity = 1;
exec = tribosolver.Execution(prec,dev,verbosity);
end

function relax = setRelaxation()
jacobianSORFactor = 6e-1;
gsSORFactor = 8e-1;
h0UpdateFactor = 5e-3;
numCycles = 3;
gamma = 2;

relax = tribosolver.Relaxation( ...
    jacobianSORFactor,gsSORFactor,h0UpdateFactor, ...
    numCycles,gamma ...
    );

end