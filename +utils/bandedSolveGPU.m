function X = bandedSolveGPU(Y,A)

% Aggregate the bands
ds = squeeze(A(2,:,:));
dl = squeeze(A(3,:,:));
d = squeeze(A(4,:,:));
du = squeeze(A(5,:,:));
dw = squeeze(A(6,:,:));
Y = squeeze(Y);

if ~isfile("+utils/bandedSolveGPUmex." + mexext)
    disp("Building cuda-mex file: Banded Solver")
    mexcuda("+utils/bandedSolveGPUmex.cu","-outdir","+utils","-lcublas")
end

ds = gpuArray(ds);
dl = gpuArray(dl);
d = gpuArray(d);
du = gpuArray(du);
dw = gpuArray(dw);
Y = gpuArray(Y);

X = utils.bandedSolveGPUmex(ds,dl,d,du,dw,Y);

X = gather(X);
end
