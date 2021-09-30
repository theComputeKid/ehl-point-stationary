function relaxPSeq( ...
    obj,l, ...
    H,hl,hr,hu,hd,hlr,hud,hudlr, ...
    dHxm3,dHxm2,dHxm1,dHx,dHxp1,dHxp2, ...
    useGS ...
    )

% Extract necessary data. Some data is read-only so there should be no
% danger of unneccesary copies.
Lk = obj.Levels(l);
nx = Lk.Domain.nx;
ny = Lk.Domain.ny;
p = Lk.Results.p;
p_rhs = Lk.p_rhs;
relGS = obj.Relaxation.gs;
relJAC = obj.Relaxation.jac;

pGS = p;
QX = zeros(nx,1,"like",p);
QY = zeros(nx,1,"like",p);
Y = zeros(nx-1,1,"like",p);
A = zeros(nx-2,6+3,"like",p);
pLine = zeros(nx,3,"like",p);

ii = 2:nx-1;
jj = 4:ny-3;
for j = jj

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

    A(1:nx-4,1) = - dHxm3(5:end,j).*iiA(5:end);

    A(1:nx-3,2) = - dHxm2(4:end,j).*iiA(4:end);
    A(1:nx-3,2) = - 0.25*hl(4:end,j).*~useGS(4:end,j).*iiA(4:end) + A(1:nx-3,2);

    A(:,3) = (hl(3:end,j) - dHxm1(3:end,j)).*iiA(3:end);
    A(:,3) = 0.25*hudlr(3:end,j).*~useGS(3:end,j).*iiA(3:end) + A(:,3);

    A(:,5) = (hr(1:end-2,j) - dHxp1(1:end-2,j)).*iiA(1:end-2);
    A(:,5) = 0.25*hudlr(1:end-2,j).*~useGS(1:end-2,j).*iiA(1:end-2) + A(:,5);

    A(4:end,6) = (- dHxp2(3:end-3,j)).*iiA(3:end-3);
    A(4:end,6) = - 0.25*hr(3:end-3,j).*~useGS(3:end-3,j).*iiA(3:end-3) + A(4:end,6);

    A(:,4) = - hudlr(ii,j) - dHx(ii,j);
    A(iiJAC-1,4) = 1.25*A(iiJAC-1,4);

    bA = flip(A.');

    kl = 3; ku = 2; transposed = false;
    b = Y(1:end-1);

    [LU, piv, ~] = matlab.internal.decomposition.builtin.bandedFactor( ...
        bA, kl, ku ...
        );

    X = [ ...
        zeros("like",p); ...
        matlab.internal.decomposition.builtin.bandedSolve(...
        LU, kl, ku, piv, b, transposed); ...
        zeros("like",p); ...
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

    if numel(pGS) < 8192
        pGS(pGS<0) = 0;
    else
        pGS(:,(j-1):(j+1)) = pGS(:,(j-1):(j+1)).*(pGS(:,(j-1):(j+1)) > 0);
    end

end

Lk.Results.p = pGS;
end