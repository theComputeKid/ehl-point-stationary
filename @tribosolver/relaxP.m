function relaxP(obj,l)

% Extract necessary data. Some data is read-only so there should be no
% danger of unneccesary copies.
Lk = obj.Levels(l);
nx = Lk.Domain.nx;
ny = Lk.Domain.ny;
dx = Lk.Domain.dx;
dy = Lk.Domain.dy;
h = Lk.Results.h;
p = Lk.Results.p;
p_rhs = Lk.p_rhs;
relGS = obj.Relaxation.gs;
relJAC = obj.Relaxation.jac;
relSwitch = obj.Relaxation.switchEpsValue;

lambda = ((128*pi^3)/(3*obj.Moes.M^4))^(1/3);
alpha = (obj.Moes.L/pi)*(3*obj.Moes.M/2)^(1/3);

lambda = lambda*exp(alpha*p);

ii = 2:nx-1;
jj = 2:ny-1;

hc = h.^3./lambda;

hd = [ ...
    zeros(1,ny,"like",h); ...
    0.5*( (h(ii + 1,:).^3)./lambda(ii + 1,:) + hc(ii,:))/dx^2; ...
    zeros(1,ny,"like",h); ...
    ];

hu = [ ...
    zeros(1,ny,"like",h); ...
    0.5*( (h(ii - 1,:).^3) ./lambda(ii - 1,:) + hc(ii,:) )./dx^2; ...
    zeros(1,ny,"like",h); ...
    ];

hr = [ ...
    zeros(nx,1,"like",h), ...
    0.5*( (h(:,jj + 1).^3)./lambda(:,jj + 1) + hc(:,jj) )./dy^2, ...
    zeros(nx,1,"like",h), ...
    ];

hl = [ ...
    zeros(nx,1,"like",h), ...
    0.5*( (h(:,jj - 1).^3) ./lambda(:,jj - 1) + hc(:,jj) )/dy^2, ...
    zeros(nx,1,"like",h) ...
    ];

H = [ ...
    zeros(1,ny,"like",h); ...
    (h(2,:) - h(1,:))/dx; ...
    (1.5*h(3:nx-1,:) - 2 * h(2:nx-2,:) + 0.5*h(1:nx-3,:))/dx; ...
    zeros(1,ny,"like",h); ...
    ];

useGS = abs(hr) > relSwitch | ...
    abs(hl) > relSwitch | ...
    abs(hu) > relSwitch | ...
    abs(hd) > relSwitch;

% Remember that we unboxed the deformation kernel
K = Lk.k(nx:nx+6,ny:ny+2);
dK0 = K(1,1) - 0.25*(2*K(2,1)+2*K(1,2)).*~useGS;
dK1 = K(2,1) - 0.25*(K(3,1)+K(1,1)+2*K(2,2)).*~useGS;
dK2 = K(3,1) - 0.25*(K(4,1)+K(2,1)+2*K(3,2)).*~useGS;
dK3 = K(4,1) - 0.25*(K(5,1)+K(3,1)+2*K(4,2)).*~useGS;
dK4 = K(5,1) - 0.25*(K(6,1)+K(4,1)+2*K(5,2)).*~useGS;

dHxm3 = (1.5*dK3 - 2*dK2 + 0.5*dK1)/dx;
dHxm2 = (1.5*dK2 - 2*dK1 + 0.5*dK0)/dx;
dHxm1 = (1.5*dK2 - 2*dK1 + 0.5*dK0)/dx;
dHx   = (1.5*dK0 - 2*dK1 + 0.5*dK2)/dx;
dHxp1 = (1.5*dK1 - 2*dK2 + 0.5*dK3)/dx;
dHxp2 = (1.5*dK2 - 2*dK3 + 0.5*dK4)/dx;

hud = hu + hd;
hlr = hl + hr;
hudlr = hud + hlr;

pGS = p;
QX = zeros(nx,1,"like",p);
QY = zeros(nx,1,"like",p);
Y = zeros(nx-1,1,"like",p);
A = zeros(nx,6,"like",p);
pLine = zeros(nx,3,"like",p);

% Calculate the number of lines that can be solved in parallel. We are
% solving a pentadiagonal system, so we need a space of 5 between lines.
%aL = 3:5:ny;

for j = jj(jj > 3 & jj < ny - 2)
    
    [iiJAC] = find(~useGS(:,j));
    [iiGS] = find(useGS(:,j));
    
    midIDx = 2;
    
    pLine(iiJAC,:) = p(iiJAC,(j-1):(j+1));
    pLine(iiGS,:) = pGS(iiGS,(j-1):(j+1));
    
    QY(ii) = ( ...
        hr(ii,j).*pLine(ii,midIDx + 1) - ...
        hlr(ii,j).*pLine(ii,midIDx) + ...
        hl(ii,j).*pLine(ii,midIDx - 1) ...
        );
    
    QX(ii) = ( ...
        hu(ii,j).*pLine(ii-1,midIDx) - ...
        hud(ii,j).*pLine(ii,midIDx) + ...
        hd(ii,j).*pLine(ii+1,midIDx) ...
        );
    
    Y(ii-1) = p_rhs(ii,j) - QX(ii) - QY(ii) + H(ii,j);
    
    iiA = [ ...
        false; ...
        pLine(ii - 1,midIDx) > 0 & pLine(ii + 1,midIDx) > 0 & ...
        pLine(ii,midIDx - 1) > 0 & pLine(ii,midIDx + 1) > 0; ...
        false; ...
        ];
    
    iiA1 = [ false; false; false; false; iiA(5:end) ];
    A(:,1) = - dHxm3(:,j).*iiA1;
    
    iiA2 = [ false; false; false; iiA(4:end) ];
    A(:,2) = - dHxm2(:,j).*iiA2;
    A(:,2) = - 0.25*hl(:,j).*~useGS(:,j).*iiA2 + A(:,2);
    
    iiA3 = [ false; false; iiA(3:end) ];
    A(:,3) = (hl(:,j) - dHxm1(:,j)).*iiA3;
    A(:,3) = 0.25*hudlr(:,j).*~useGS(:,j).*iiA3 + A(:,3);
    
    iiA5 = [ iiA(1:end-1); false ];
    A(:,5) = (hr(:,j) - dHxp1(:,j)).*iiA5;
    A(:,5) = 0.25*hudlr(:,j).*~useGS(:,j).*iiA5 + A(:,5);
    
    iiA6 = [ iiA(1:end-1); false ];
    A(:,6) = (- dHxp2(:,j)).*iiA6;
    A(:,6) = - 0.25*hr(:,j).*~useGS(:,j).*iiA6 + A(:,6);
    
    A(ii,4) = - hudlr(ii,j) - dHx(ii,j);
    A(iiJAC,4) = 1.25*A(iiJAC,4);
    
    AA = zeros(3+size(A,2),size(A,1)-2,"like",A);
    AA(4,3:end) = A(2:nx-3,6);
    AA(5,2:end) = A(2:nx-2,5);
    AA(6,:) = A(2:nx-1,4);
    AA(7,1:end-1) = A(3:nx-1,3);
    AA(8,1:end-2) = A(4:nx-1,2);
    AA(9,1:end-3) = A(5:nx-1,1);
    
    kl = 3;
    ku = 2;
    transposed = false;
    b = Y(1:end-1);
    [LU, piv, ~] = matlab.internal.decomposition.builtin.bandedFactor( ...
        AA, kl, ku ...
        );
    
    X = [ ...
        0; ...
        matlab.internal.decomposition.builtin.bandedSolve(...
        LU, kl, ku, piv, b, transposed); ...
        0; ...
        ];
    
    p0 = pGS(:,j);
    pGS(iiGS,j) = pGS(iiGS,j) + X(iiGS)*relGS;
    pGS(iiJAC,j) = pGS(iiJAC,j) + X(iiJAC)*relJAC;
    del = pGS(:,j) - p0;
    
    iiJAC(iiJAC < 3 | iiJAC > (nx - 2)) = [];
    pGS(iiJAC - 1, j) = pGS(iiJAC - 1, j) - 0.25*del(iiJAC);
    pGS(iiJAC, j - 1) = pGS(iiJAC, j - 1) - 0.25*del(iiJAC);
    pGS(iiJAC + 1, j) = pGS(iiJAC + 1, j) - 0.25*del(iiJAC);
    pGS(iiJAC, j + 1) = pGS(iiJAC, j + 1) - 0.25*del(iiJAC);
    
    pGS(pGS<0) = 0;
    
end

Lk.Results.p = pGS;

Lk.calcDeformation();

if obj.Execution.Verbosity > 1
    surf(Lk.Domain.x,Lk.Domain.y,pGS);
    xlabel("X"); ylabel("Y"); zlabel("P");
    title("P: Level " + l)
    drawnow
end
end