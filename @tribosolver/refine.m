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

iic = 2:nxc;
iif = iif(2:end);
jjc = 2:nyc;
jjf = jjf(2:end);

p(iif-1,jjf) = p(iif-1,jjf) + pcim1(iic,jjc).*(p(iif-1,jjf) > 0);
p(iif,jjf-1) = p(iif,jjf-1) + pcjm1(iic,jjc).*(p(iif,jjf-1) > 0);
p(iif-1,jjf-1) = p(iif-1,jjf-1) + pcim1jm1(iic,jjc).*(p(iif-1,jjf-1) > 0);

p(p<0)=0;
obj.Levels(l).Results.p = p;
end
