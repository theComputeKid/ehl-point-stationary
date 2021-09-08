function refine(obj,l)

Lk = obj.Levels(l);
Lkc = obj.Levels(l-1);

nx = Lk.Domain.nx;
nxc = Lkc.Domain.nx;

ny = Lk.Domain.ny;
nyc = Lkc.Domain.ny;

p = Lk.Results.p;
pc = Lkc.Results.p;
pco = Lkc.p_old;

iif = 1:2:nx;
jjf = 1:2:ny;

p(iif,jjf) = p(iif,jjf) + (pc - pco).*(p(iif,jjf) > 0);

pcim1(2:nxc,:) = ( ...
    pc(2:end,:) - pco(2:end,:) + ...
    pc(1:end-1,:) - pco(1:end-1,:) ...
    )*0.5;

pcjm1(:,2:nyc) = ( ...
    pc(:,2:nyc) - pco(:,2:nyc) + ...
    pc(:,1:end-1) - pco(:,1:end-1) ...
    )*0.5;

pcim1jm1(2:nxc,2:nyc) = ( ...
    pc(2:nxc,2:nyc) - pco(2:nxc,2:nyc) + ...
    pc(2:nxc,1:end-1) - pco(2:nxc,1:end-1) + ...
    pc(1:end-1,2:nyc) - pco(1:end-1,2:nyc) + ...
    pc(1:end-1,1:end-1) - pco(1:end-1,1:end-1) ...
    )*0.25;

[pcim1, pcjm1, pcim1jm1, p] = gather(pcim1,pcjm1,pcim1jm1, p);

for i=2:nxc
    
    iff = 2*i - 1;
    
    for j=2:nyc
        
        jff = 2*j - 1;
        
        if p(iff - 1,jff) > 0 && j < nyc
            p(iff - 1,jff) = p(iff - 1,jff) + ...
                pcim1(i,j);
        end
        
        if p(iff, jff - 1) > 0 && i < nxc
            p(iff, jff - 1) = p(iff, jff - 1) + ...
                pcjm1(i,j);
        end
        
        if (p(iff - 1, jff - 1) > 0)
            p(iff - 1, jff - 1) = p(iff - 1, jff - 1) + ...
                pcim1jm1(i,j);
        end
    end
end

if isgpuarray(pc)
    p = gpuArray(p);
end

p(p<0)=0;
obj.Levels(l).Results.p = p;
end
