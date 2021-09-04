function results = solve(obj)
% SOLVE:
% Main solution driver routine for tribosolver

obj.initLevels();
obj.mgFull(obj.Domain.mgl);
results = obj.Levels(obj.Domain.mgl).Results;
end