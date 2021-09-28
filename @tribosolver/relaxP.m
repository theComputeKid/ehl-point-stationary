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

obj.relaxPSeq( ...
    l, ...
    H,hl,hr,hu,hd,hlr,hud,hudlr, ...
    dHxm3,dHxm2,dHxm1,dHx,dHxp1,dHxp2, ...
    useGS ...
    );

Lk.calcDeformation();

if obj.Execution.Verbosity > 1
    surf(Lk.Domain.x,Lk.Domain.y,pGS);
    xlabel("X"); ylabel("Y"); zlabel("P");
    title("P: Level " + l)
    drawnow
end
end