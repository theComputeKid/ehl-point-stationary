function results = solve(obj)
% SOLVE:
% Main solution driver routine for tribosolver

obj.initLevels();

for i = 1:obj.Relaxation.numCycles
    obj.mgFull(obj.Domain.mgl);
end

results = obj.Levels(obj.Domain.mgl).Results;
end