function mgCycle(obj,l)
% Executes one complete V or W multigrid cycle

if (l == 1)
    
    for i = 1:obj.Relaxation.itMain
        obj.relax(l);
    end
    
else
    
    for i = 1:obj.Relaxation.itPre
        obj.relax(l);
    end
    
    obj.coarsen(l);
    
    for i=1:obj.Relaxation.gamma
        obj.mgCycle(l - 1);
    end
    
    obj.refine(l);
    
    for i = 1:obj.Relaxation.itPost
        obj.relax(l);
    end
    
end

end