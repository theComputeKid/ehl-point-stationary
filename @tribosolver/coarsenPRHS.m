function coarsenPRHS(obj,l)

Lk = obj.Levels(l);
Lkc = obj.Levels(l-1);

dx = Lk.Domain.dx;
dxc = Lkc.Domain.dx;

dy = Lk.Domain.dy;
dyc = Lkc.Domain.dy;

nxc = Lkc.Domain.nx;
nyc = Lkc.Domain.ny;

p = Lk.Results.p;
pc = Lkc.Results.p;

p_rhs = Lk.p_rhs;
p_rhsc = Lkc.p_rhs;
h = Lk.Results.h;
hc = Lkc.Results.h;

iic=2:nxc-1; jjc=2:nyc-1;
iif=2*iic-1; jjf=2*jjc-1;

lambda = ((128*pi^3)/(3*obj.Moes.M^4))^(1/3);

LUf = LU(h,p,dx,dy,lambda);
c(iif,jjf) = p_rhs(iif,jjf)-LUf(iif,jjf);
n(iic,jjc) = (p_rhs(iif,jjf + 1) - LUf(iif, jjf + 1)).*(p(iif,jjf + 1) > 0);
e(iic,jjc) = (p_rhs(iif + 1,jjf) - LUf(iif + 1, jjf)).*(p(iif + 1,jjf) > 0);
s(iic,jjc) = (p_rhs(iif,jjf - 1) - LUf(iif, jjf - 1)).*(p(iif,jjf - 1) > 0);
w(iic,jjc) = (p_rhs(iif - 1,jjf) - LUf(iif - 1, jjf)).*(p(iif - 1,jjf) > 0);
ne(iic,jjc) = (p_rhs(iif + 1,jjf + 1) - LUf(iif + 1, jjf + 1)).*(p(iif + 1,jjf + 1) > 0);
nw(iic,jjc) = (p_rhs(iif - 1,jjf + 1) - LUf(iif - 1, jjf + 1)).*(p(iif - 1,jjf + 1) > 0);
se(iic,jjc) = (p_rhs(iif + 1,jjf - 1) - LUf(iif + 1, jjf - 1)).*(p(iif + 1,jjf - 1) > 0);
sw(iic,jjc) = (p_rhs(iif - 1,jjf - 1) - LUf(iif - 1, jjf - 1)).*(p(iif - 1,jjf - 1) > 0);

p_rhsc(iic,jjc) = (4 * c(iif,jjf) + 2 * (n(iic,jjc) + s(iic,jjc) + e(iic,jjc) + w(iic,jjc)) + (ne(iic,jjc) + se(iic,jjc) + sw(iic,jjc) + nw(iic,jjc))) / 16.0;
p_rhsc = p_rhsc + LU(hc,pc,dxc,dyc,lambda);

obj.Levels(l-1).p_rhs = p_rhsc;

end



%% Linear Operator
function out = LU(h,p,dx,dy,lambda)
nx = size(h,1); ny = size(h,2); ii = 2:nx-1; jj = 2:ny - 1;

hc = zeros(nx,ny);
hu = hc; hd = hc; hl = hc; hr = hc; H = hc; QX=H; QY=H;

hc(ii,jj) = h(ii,jj).^3;

hr(ii,jj) = 0.5*(h(3:nx,jj).^3 + hc(ii,jj))/lambda./dx.^2;
hl(ii,jj) = 0.5*(h(1:nx-2,jj).^3 + hc(ii,jj))/lambda./dx.^2;

hu(ii,jj) = 0.5*(h(ii,jj + 1).^3 + hc(ii,jj))/lambda./dy.^2;
hd(ii,jj) = 0.5*(h(ii,jj - 1).^3 + hc(ii,jj))/lambda./dy.^2;

H(2,jj) = (h(2,jj) - h(1,jj))./dx;
H(3:nx-1,jj) = (1.5*h(3:nx-1,jj) - 2 * h(2:nx-2,jj) + 0.5*h(1:nx-3,jj))./dx;

QX(ii,jj) = hr(ii,jj).*p(ii+1,jj) - (hr(ii,jj) + hl(ii,jj)).*p(ii,jj) + hl(ii,jj).*p(ii-1,jj);
QY(ii,jj) = hu(ii,jj).*p(ii,jj+1) - (hu(ii,jj) + hd(ii,jj)).*p(ii,jj) + hd(ii,jj).*p(ii,jj-1);

out = (QX+QY-H);
end