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

    iiA = cat( 1, ...
        false(1,1,cnL), ...
        pLine(ii - 1,midIDx,:) > 0 & pLine(ii + 1,midIDx,:) > 0 & ...
        pLine(ii,midIDx - 1,:) > 0 & pLine(ii,midIDx + 1,:) > 0, ...
        false(1,1,cnL) ...
        );

    iiA1 = cat( 1, ...
        false(4,1,cnL), ...
        iiA(5:end,1,:) ...
        );

    A(:,1,1:cnL) = reshape(-dHxm3(:,caL),[],1,cnL).*iiA1;

    iiA2 = cat( 1, ...
        false(3,1,cnL), ...
        iiA(4:end,1,:) ...
        );

    A(:,2,1:cnL) = reshape(-dHxm2(:,caL),[],1,cnL).*iiA2;
    A(:,2,1:cnL) = reshape(-0.25*hl(:,caL),[],1,cnL).*~aGS.*iiA2 + A(:,2,1:cnL);

    iiA3 = cat( 1, ...
        false(2,1,cnL), ...
        iiA(3:end,1,:) ...
        );

    A(:,3,1:cnL) = reshape((hl(:,caL) - dHxm1(:,caL)),[],1,cnL).*iiA3;
    A(:,3,1:cnL) = reshape(0.25*hudlr(:,caL),[],1,cnL).*~aGS.*iiA3 + A(:,3,1:cnL);

    iiA5 = cat( 1, ...
        iiA(1:end-1,1,:), ...
        false(1,1,cnL) ...
        );

    A(:,5,1:cnL) = reshape((hr(:,caL) - dHxp1(:,caL)),[],1,cnL).*iiA5;
    A(:,5,1:cnL) = reshape(0.25*hudlr(:,caL),[],1,cnL).*~aGS.*iiA5 + A(:,5,1:cnL);

    iiA6 = iiA5;
    A(:,6,1:cnL) = reshape((-dHxp2(:,caL)),[],1,cnL).*iiA6;
    A(:,6,1:cnL) = reshape(-0.25*hr(:,caL),[],1,cnL).*~aGS.*iiA6 + A(:,6,1:cnL);

    A(ii,4,1:cnL) = - hudlr(ii,caL) - dHx(ii,caL);
    A(~aGS) = 1.25*A(~aGS);

    for n = 1:cnL
        AA(4,3:end,n) = A(2:nx-3,6,n);
        AA(5,2:end,n) = A(2:nx-2,5,n);
        AA(6,:,n) = A(2:nx-1,4,n);
        AA(7,1:end-1,n) = A(3:nx-1,3,n);
        AA(8,1:end-2,n) = A(4:nx-1,2,n);
        AA(9,1:end-3,n) = A(5:nx-1,1,n);
    end

    if obj.Execution.Device == "cpu_par"
        X = bandedSolverCPU(Y,AA,cnL);
    elseif obj.Execution.Device == "gpu"
        X = utils.bandedSolveGPU(Y,AA(:,:,1:cnL));
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

function X = bandedSolverCPU(Y,A,cnL)

nx = size(A,2) + 2; 
X = zeros(nx,1,cnL,"like",Y);
parfor n = 1:cnL
    kl = 3; ku = 2; transposed = false;
    b = Y(:,1,n);
    [LU, piv, ~] = matlab.internal.decomposition.builtin.bandedFactor( ...
        A(:,:,n), kl, ku ...
        );

    X(:,:,n) = [ ...
        0; ...
        matlab.internal.decomposition.builtin.bandedSolve(...
        LU, kl, ku, piv, b(1:end-1), transposed); ...
        0; ...
        ];
end

end
