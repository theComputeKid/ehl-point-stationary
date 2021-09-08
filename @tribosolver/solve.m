function results = solve(obj)
% SOLVE:
% Main solution driver routine for tribosolver

obj.initLevels();

if obj.Execution.Verbosity
    disp("Starting solution scheme. Convergence monitor:")
end

tic
for i = 1:obj.Relaxation.numCycles
    p = obj.Levels(obj.Domain.mgl).Results.p;
    h = obj.Levels(obj.Domain.mgl).Results.h;
    obj.mgFull(obj.Domain.mgl);
    if checkTol(obj,p,h,i)
        break;
    end
end

t = toc;
if obj.Execution.Verbosity
    disp("Time to solve (s): " + t)
end

results = obj.Levels(obj.Domain.mgl).Results;
end

function hasConverged = checkTol(obj,p,h,N)
pNew = obj.Levels(obj.Domain.mgl).Results.p;
hNew = obj.Levels(obj.Domain.mgl).Results.h;

pTol = max(abs(pNew - p),[],"all");
hTol = max(abs(hNew - h),[],"all");

fbNew = sum( ...
    pNew .* ...
    obj.Levels(obj.Domain.mgl).Domain.dx .* ...
    obj.Levels(obj.Domain.mgl).Domain.dy, ...
    "all");

fb = sum( ...
    p .* ...
    obj.Levels(obj.Domain.mgl).Domain.dx .* ...
    obj.Levels(obj.Domain.mgl).Domain.dy, ...
    "all");

fbTol = abs(fbNew - fb);

if obj.Execution.Verbosity
    disp("N: " + N + " H0: " + obj.Levels(obj.Domain.mgl).Results.h0 + ...
        " dP: " + pTol + " dH: " + hTol + " dFb: " + fbTol)
end

hasConverged = (pTol < obj.Tolerance.p) && ...
    (hTol < obj.Tolerance.h) && ...
    (fbTol < obj.Tolerance.fb);

end