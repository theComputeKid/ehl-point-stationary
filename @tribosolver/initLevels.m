function initLevels(obj)
% Initialize the multigrid levels, allocating space and setting initial
% values.

mgl = obj.Domain.mgl;

% Pre-allocate structs
obj.Levels(1:mgl) = tribosolver.internal.Level;

for k = 1:mgl
    obj.Levels(k) = tribosolver.internal.Level( ...
        obj.Domain, obj.Execution, k ...
        );
end

end