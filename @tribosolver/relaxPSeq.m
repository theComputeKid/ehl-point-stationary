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
A = zeros(nx,6,"like",p);
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
end