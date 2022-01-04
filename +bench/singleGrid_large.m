function results = singleGrid_large(backend,precision)

% singleGrid_large: Benchmarks Single Grid solver. Reasonably long running.
%
% To run, go to the project root directory and type: bench.singleGrid_large
% Optionally specify backend (cpu_seq, cpu_par, gpu_par) and precision
% (single, double). E.g.: bench.singleGrid_large("cpu_par","single").
% Defaults to "cpu_seq" and "double".
%
% Copyright (C) 2021-2022 theComputeKid

narginchk(0,2)

if ~nargin
    backend = "cpu_seq";
end

if nargin < 2
    precision = "double";
end

domain = setDomain();
moes = setMoes();
exec = setExecution(backend,precision);
relax = setRelaxation();

% Ensure parallel pool creation not timed as part of benchmarks.
if backend == "cpu_par"
    if isempty(gcp('nocreate'))
        parpool("threads");
    end
end

disp("Backend: " + backend)
disp("Precision: " + precision)

if backend == "cpu_par" || backend == "cpu_seq"
    info = utils.cpuinfo;
    disp(info.CPUName)
else
    disp(gpuDevice().Name)
end

model = tribosolver(relax,domain,exec,moes);
results = model.solve();

end

function domain = setDomain()

nx = 4096; xin = -2.5; xout = 1.5;
ny = 4096; yin = -2.5; yout = 2.5;
mgl = 1;

domain = tribosolver.Domain(xin,xout,nx,yin,yout,ny,mgl);

end

function moes = setMoes()
M = 15; L = 0;
H0 = -0.63;
moes = tribosolver.Moes(M,L,H0);
end

function exec = setExecution(backend,precision)
verbosity = 1;
exec = tribosolver.Execution(precision,backend,verbosity);
end

function relax = setRelaxation()

jacobianSORFactor = 6e-1;
gsSORFactor = 8e-1;
h0UpdateFactor = 1e-3;
epsSwitch = 0.3;

numCycles = 10;
gamma = 2;

itPre = 3;
itPost = 3;
itMain = 5;

relax = tribosolver.Relaxation( ...
    jacobianSORFactor,gsSORFactor,h0UpdateFactor, ...
    numCycles,gamma, ...
    epsSwitch, ...
    itPre, itMain, itPost ...
    );

end
