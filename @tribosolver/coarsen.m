function coarsen(obj,l)

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

for i=1:nxc
    iff = 2*i - 1;
    for j=1:nyc
        jff = 2*j - 1;
        
        if i==1 || j==1 || i == nxc || j==nyc
            pc(i,j)=0;
        else
            if ((p(iff ,jff ) == 0) || ...
                (p(iff + 1,jff + 1) == 0) || (p(iff + 1,jff - 1) == 0) || ...
                (p(iff - 1,jff + 1) == 0) || (p(iff - 1,jff - 1) == 0) || ...
                (p(iff + 1,jff ) == 0) || (p(iff - 1,jff ) == 0) || ...
                (p(iff ,jff + 1) == 0) || (p(iff ,jff - 1) == 0))

                pc(i,j) = p(iff,jff); % Injection if near cavitation point       
            else
                pc(i,j) =(4.0*p(iff,jff) + ...
                2.0*(p(iff + 1,jff) + p(iff - 1,jff) + p(iff,jff + 1) + p(iff,jff - 1)) + ...
                p(iff + 1,jff + 1) + p(iff + 1,jff - 1) + p(iff - 1,jff + 1) + p(iff - 1,jff - 1)) / 16.0;
            end
        end
    end
end

obj.Levels(l-1).Results.p = pc;
obj.Levels(l-1).p_old = pc;

%% Coarsen P RHS
p_rhs = Lk.p_rhs;
p_rhsc = Lkc.p_rhs;
h = Lk.Results.h;
hc = Lkc.Results.h;

iic=2:nxc-1; jjc=2:nyc-1;
iif=2*iic-1; jjf=2*jjc-1;

lambda = ((128*pi^3)/(3*obj.Moes.M^4))^(1/3);

LUf = LU(h,p,dx,dy,lambda);
c(iif,jjf) = p_rhs(iif,jjf)-LUf(iif,jjf);

for i=iic
    iff = 2*i - 1;
    for j=jjc
        jff = 2*j - 1;
                
        if (p(iff,jff + 1) > 0)
            n = p_rhs(iff,jff + 1) - LUf(iff, jff + 1);
        else
            n = 0;
        end
        
        if (p(iff + 1,jff) > 0)
            e = p_rhs(iff + 1,jff) - LUf(iff + 1, jff);
        else
            e = 0;
        end
        
        if (p(iff,jff - 1) > 0)
            s = p_rhs(iff,jff - 1) - LUf(iff, jff - 1);
        else
            s = 0;
        end
        
        if (p(iff - 1,jff) > 0)
            w = p_rhs(iff - 1,jff) - LUf(iff - 1, jff);
        else
            w = 0;
        end
        
        if (p(iff + 1,jff + 1) > 0)
            ne = p_rhs(iff + 1,jff + 1) - LUf(iff + 1, jff + 1);
        else
            ne = 0;
        end
        
        if (p(iff - 1,jff + 1) > 0)
            nw = p_rhs(iff - 1,jff + 1) - LUf(iff - 1, jff + 1);
        else
            nw = 0;
        end
        
        if (p(iff + 1,jff - 1) > 0)
            se = p_rhs(iff + 1,jff - 1) - LUf(iff + 1, jff - 1);
        else
            se = 0.0;
        end
        
        if (p(iff - 1,jff - 1) > 0)
            sw = p_rhs(iff - 1,jff - 1) - LUf(iff - 1, jff - 1);
        else
            sw = 0.0;
        end

        p_rhsc(i,j) = (4 * c(iff,jff) + 2 * (n + s + e + w) + (ne + se + sw + nw)) / 16.0;
    end
end

p_rhsc = p_rhsc + LU(hc,pc,dxc,dyc,lambda);
obj.Levels(l-1).p_rhs = p_rhsc;

%% Coarsen H
obj.Levels(l-1).Results.h = downsample(downsample(h,2)',2)';

%% Coarsen Force Balance
Ff = sum(p.*dx.*dy,"all");
Fc = sum(pc.*dxc.*dyc,"all");

obj.Levels(l-1).fb = Lk.fb + Ff - Fc;
end

%% Linear Operator
function out = LU(h,p,dx,dy,lambda)
    nx = size(h,1); ny = size(h,2); ii = 2:nx-1; jj = 2:ny - 1;
    
    hc = zeros(nx,ny);    
    hu = hc; hd = hc; hl = hc; hr = hc; H = hc; QX=H; QY=H;
    
    hc(ii,jj) = h(ii,jj).^3;

    hr(ii,jj) = 0.5*(h(3:nx,jj).^3 + hc(ii,jj))/lambda/dx^2;
    hl(ii,jj) = 0.5*(h(1:nx-2,jj).^3 + hc(ii,jj))/lambda/dx^2; 

    hu(ii,jj) = 0.5*(h(ii,jj + 1).^3 + hc(ii,jj))/lambda/dy^2; 
    hd(ii,jj) = 0.5*(h(ii,jj - 1).^3 + hc(ii,jj))/lambda/dy^2;
    
    H(2,jj) = (h(2,jj) - h(1,jj))/dx;
    H(3:nx-1,jj) = (1.5*h(3:nx-1,jj) - 2 * h(2:nx-2,jj) + 0.5*h(1:nx-3,jj))/dx;
    
    QX(ii,jj) = hr(ii,jj).*p(ii+1,jj) - (hr(ii,jj) + hl(ii,jj)).*p(ii,jj) + hl(ii,jj).*p(ii-1,jj);
    QY(ii,jj) = hu(ii,jj).*p(ii,jj+1) - (hu(ii,jj) + hd(ii,jj)).*p(ii,jj) + hd(ii,jj).*p(ii,jj-1);
        
    out = (QX+QY-H);
end