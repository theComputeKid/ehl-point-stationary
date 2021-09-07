function mgFull(obj,l)
% Driver for the full multigrid cycle starting at level l

if (l == 1)
    
    if obj.Domain.mgl == 1
        
        for i = 1:obj.Relaxation.numCycles
            obj.relaxP(l);
            obj.relaxH0(l);
        end
        
    else
        
        for i = 1:obj.Relaxation.itMain
            obj.relaxP(l);
            obj.relaxH0(l);
        end
        
    end
    
else
    
    obj.mgFull(l - 1);
    
    obj.mgFullInterp(l);
    
    for i=1:obj.Relaxation.numCycles
        obj.mgCycle(l);
    end
    
end

end