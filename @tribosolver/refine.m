function refine(obj,l)

Lk = obj.Levels(l);
Lkc = obj.Levels(l-1);

dx = Lk.Domain.dx;
dxc = Lkc.Domain.dx;

dy = Lk.Domain.dy;
dyc = Lkc.Domain.dy;

nx = Lk.Domain.nx;
nxc = Lkc.Domain.nx;

ny = Lk.Domain.ny;
nyc = Lkc.Domain.ny;

p = Lk.Results.p;
pc = Lkc.Results.p;
pco = Lkc.p_old;

%% Refine P
for i=2:nxc
    iff = 2*i - 1;
    for j=2:nyc
        jff = 2*j - 1;
       
        if p(iff,jff) > 0
            p(iff,jff) = p(iff,jff) + (pc(i,j) - pco(i,j));
        end
        
        if p(iff - 1,jff) > 0 && j < nyc
            p(iff - 1,jff) = p(iff - 1,jff) + ...
                (pc(i,j) - pco(i,j) + pc(i - 1,j) - pco(i - 1,j))*0.5;
        end
        
        if p(iff,jff - 1)>0 && i < nxc
            p(iff,jff - 1) = p(iff,jff - 1) + ...
                (pc(i,j) - pco(i,j) + pc(i,j - 1) - pco(i,j - 1))*0.5;
        end
        
        if (p(iff - 1,jff - 1)>0)
            p(iff - 1,jff - 1) = p(iff - 1,jff - 1)+ ...
                (pc(i,j) - pco(i,j) + pc(i,j - 1) - pco(i,j - 1) + ...
            pc(i - 1,j) - pco(i - 1,j) + pc(i - 1,j - 1) - pco(i - 1,j - 1))*0.25;
        end
    end
end

p(p<0)=0;
obj.Levels(l).Results.p = p;

end