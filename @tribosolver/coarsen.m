function coarsen(obj,l)

coarsenP(obj,l);
coarsenPRHS(obj,l);
coarsenFB(obj,l);

obj.Levels(l-1).Results.h = downsample( ...
    downsample(obj.Levels(l).Results.h,2)' ...
    ,2)';

end

