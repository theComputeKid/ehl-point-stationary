function relaxPPar( ...
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

% Calculate the number of lines that can be solved in parallel. We are
% solving a sexadiagonal system, so we need a space of 6 between lines.
nD = 6; % Number of diagonals being solved
ii = 2:nx-1;
jj = 4:ny-3;
aL = jj(1):nD:jj(end); % Currently active lines being solved
nL = length(aL); % Number of active lines

Y = zeros(nx-1,1,nL,"like",p);
A = zeros(nx,6,nL,"like",p);
AA = zeros(3+size(A,2),size(A,1)-2,nL,"like",A);

for it = 0:(nD-1)

    % Currently active lines. Ensure that they are not out-of-bounds
    caL = aL + it;
    caL = caL(caL < aL(end));

    % Whether to use GS/JAC for each line
    aGS = useGS(:,caL);
    aGS = reshape(aGS,nx,1,[]);
    cnL = length(caL);

    aLRegion = reshape(reshape([caL-1, caL, caL+1],cnL,3)',1,[]);
    repGS = reshape(repelem(aGS,3,1),[],cnL);

    pLine = reshape( ...
        reshape(p(:,aLRegion),[],cnL).*repGS + ...
        reshape(pGS(:,aLRegion),[],cnL).*~repGS, ...
        [],3,cnL);

    midIDx = 2;

    QY = ( ...
        reshape(hr(:,caL),[],1,cnL).*pLine(:,midIDx + 1,:) - ...
        reshape(hlr(:,caL),[],1,cnL).*pLine(:,midIDx,:) + ...
        reshape(hl(:,caL),[],1,cnL).*pLine(:,midIDx - 1,:) ...
        );

    QX(ii,1,1:cnL) = ( ...
        reshape(hu(ii,caL),[],1,cnL).*pLine(ii-1,midIDx,:) - ...
        reshape(hud(ii,caL),[],1,cnL).*pLine(ii,midIDx,:) + ...
        reshape(hd(ii,caL),[],1,cnL).*pLine(ii+1,midIDx,:) ...
        ); %#ok<AGROW>

    Y(ii-1,1,1:cnL) = reshape(p_rhs(ii,caL),[],1,cnL) - ...
        QX(ii,1,1:cnL) - QY(ii,1,1:cnL) + ...
        reshape(H(ii,caL),[],1,cnL);

    for n = 1:cnL
        ll = caL(n);

        iiA = [ ...
            false; ...
            pLine(ii - 1,midIDx,n) > 0 & pLine(ii + 1,midIDx,n) > 0 & ...
            pLine(ii,midIDx - 1,n) > 0 & pLine(ii,midIDx + 1,n) > 0; ...
            false; ...
            ];

        iiA1 = [ false; false; false; false; iiA(5:end) ];
        A(:,1,n) = - dHxm3(:,ll).*iiA1;

        iiA2 = [ false; false; false; iiA(4:end) ];
        A(:,2,n) = - dHxm2(:,ll).*iiA2;
        A(:,2,n) = - 0.25*hl(:,ll).*~useGS(:,ll).*iiA2 + A(:,2,n);

        iiA3 = [ false; false; iiA(3:end) ];
        A(:,3,n) = (hl(:,ll) - dHxm1(:,ll)).*iiA3;
        A(:,3,n) = 0.25*hudlr(:,ll).*~useGS(:,ll).*iiA3 + A(:,3,n);

        iiA5 = [ iiA(1:end-1); false ];
        A(:,5,n) = (hr(:,ll) - dHxp1(:,ll)).*iiA5;
        A(:,5,n) = 0.25*hudlr(:,ll).*~useGS(:,ll).*iiA5 + A(:,5,n);

        iiA6 = [ iiA(1:end-1); false ];
        A(:,6,n) = (- dHxp2(:,ll)).*iiA6;
        A(:,6,n) = - 0.25*hr(:,ll).*~useGS(:,ll).*iiA6 + A(:,6,n);

        A(ii,4,n) = - hudlr(ii,ll) - dHx(ii,ll);
        A(~useGS(:,ll),4,n) = 1.25*A(~useGS(:,ll),4,n);

        AA(4,3:end,n) = A(2:nx-3,6,n);
        AA(5,2:end,n) = A(2:nx-2,5,n);
        AA(6,:,n) = A(2:nx-1,4,n);
        AA(7,1:end-1,n) = A(3:nx-1,3,n);
        AA(8,1:end-2,n) = A(4:nx-1,2,n);
        AA(9,1:end-3,n) = A(5:nx-1,1,n);
    end

    X = zeros(nx,1,cnL);
    parfor n = 1:cnL
        kl = 3; ku = 2; transposed = false;
        b = Y(:,1,n);
        [LU, piv, ~] = matlab.internal.decomposition.builtin.bandedFactor( ...
            AA(:,:,n), kl, ku ...
            );

        X(:,:,n) = [ ...
            0; ...
            matlab.internal.decomposition.builtin.bandedSolve(...
            LU, kl, ku, piv, b(1:end-1), transposed); ...
            0; ...
            ];
    end

    p0 = pGS(:,caL);


    pGS(:,caL) = pGS(:,caL) + ...
        reshape( ...
        X.*(relGS.*aGS+relJAC.*~aGS), ...
        nx,[]);
    del = pGS(:,caL) - p0;

    for n = 1:cnL
        ll = caL(n);
        [iiJAC] = find(~useGS(:,ll));

        iiJAC(iiJAC < 3 | iiJAC > (nx - 2)) = [];
        pGS(iiJAC - 1, ll) = pGS(iiJAC - 1, ll) - 0.25*del(iiJAC,n);
        pGS(iiJAC, ll - 1) = pGS(iiJAC, ll - 1) - 0.25*del(iiJAC,n);
        pGS(iiJAC + 1, ll) = pGS(iiJAC + 1, ll) - 0.25*del(iiJAC,n);
        pGS(iiJAC, ll + 1) = pGS(iiJAC, ll + 1) - 0.25*del(iiJAC,n);

    end

    pGS(pGS<0) = 0;

end
Lk.Results.p = pGS;
end