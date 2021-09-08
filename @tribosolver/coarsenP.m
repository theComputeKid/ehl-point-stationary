function coarsenP(obj,l)

Lk = obj.Levels(l);
Lkc = obj.Levels(l-1);

nx = Lk.Domain.nx;
nxc = Lkc.Domain.nx;

ny = Lk.Domain.ny;
nyc = Lkc.Domain.ny;

p = Lk.Results.p;
pc = Lkc.Results.p;

iif = 3:2:nx - 1;
jjf = 3:2:ny - 1;

mask = false(size(p));

mask(iif,jjf) = (p(iif, jjf) == 0) | ...
    (p(iif + 1, jjf + 1) == 0) | (p(iif + 1,jjf - 1) == 0) | ...
    (p(iif - 1, jjf + 1) == 0) | (p(iif - 1,jjf - 1) == 0) | ...
    (p(iif + 1, jjf) == 0) | (p(iif - 1,jjf ) == 0) | ...
    (p(iif, jjf + 1) == 0) | (p(iif ,jjf - 1) == 0);

pctemp = zeros(nxc,nyc,"like",pc);
pctemp(2:nxc-1,2:nyc-1) = (4.0*p(iif,jjf) + ...
    2.0*(p(iif + 1,jjf) + p(iif - 1,jjf) + p(iif,jjf + 1) + p(iif,jjf - 1)) + ...
    p(iif + 1,jjf + 1) + p(iif + 1,jjf - 1) + p(iif - 1,jjf + 1) + p(iif - 1,jjf - 1)) / 16.0;

pc = p(1:2:nx,1:2:ny).*mask(1:2:nx,1:2:ny) + ...
    pctemp.*~mask(1:2:nx,1:2:ny);

obj.Levels(l-1).Results.p = pc;
obj.Levels(l-1).p_old = pc;

end