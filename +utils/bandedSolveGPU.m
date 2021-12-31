function X = bandedSolveGPU(A,Y)

% Aggregate the bands
ds = squeeze(A(:,2,:));
dl = squeeze(A(:,3,:));
d =  squeeze(A(:,4,:));
du = squeeze(A(:,5,:));
dw = squeeze(A(:,6,:));
Y = Y(1:end-1,:);

if ~isfile("+utils/pentasolver." + mexext)
    disp("Building cuda-mex file: Banded Solver")
    if ispc
        !cd +utils && .\make-mex.bat
    end
end

X = zeros(size(Y,1) + 2, size(Y,2),"like",Y);

X(2:end-1,:) = utils.pentasolver(Y,ds,dl,d,du,dw);

X = reshape(X,size(X,1),1,size(X,2));
end
