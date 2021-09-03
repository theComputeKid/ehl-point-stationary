function results = solve(obj)
% SOLVE:
% Main solution driver routine for tribosolver

obj.initLevels();

for kc = 1:obj.Relaxation.numCycles
    
    obj.mgCycle();
    
end

results = tribosolver.Results;
end