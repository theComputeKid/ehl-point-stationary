function mgFull(obj,l)
% Driver for the full multigrid cycle starting at level l

if (l == 1)
    
    for i = 1:obj.Relaxation.itMain
        obj.relaxP(l);
        obj.relaxH0(l);
    end
    
else
    
    obj.mgFull(l - 1);
    
    obj.mgFullInterp(l);
    
    for i=1:obj.Relaxation.gamma
        obj.mgCycle(l);
    end
    
end

end