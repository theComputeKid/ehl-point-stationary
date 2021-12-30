function X = bandedSolveGPU(A,Y)

% Aggregate the bands
ds = A(:,2);
dl = A(:,3);
d =  A(:,4);
du = A(:,5);
dw = A(:,6);
Y = Y(1:end-1,:);

if ~isfile("+utils/pentasolver." + mexext)
    disp("Building cuda-mex file: Banded Solver")
    if ispc
        !cd +utils && nmake
    end
end

X = [ ...
    zeros("like",Y); ...
    utils.pentasolver(Y,ds,dl,d,du,dw); ...
    zeros("like",Y); ...
    ];

end
