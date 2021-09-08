classdef Level < handle
    
    % Container for a single MultiGrid Level. Internal solver use only.
    %
    % Copyright (C) 2021 theComputeKid
    
    properties(SetAccess=immutable)
        
        Domain(1,1) tribosolver.internal.Domain
        
        k(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
        k_fft(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
        h(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
    end
    
    properties
        
        Results(1,1) tribosolver.Results
        
        p_rhs(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
        p_old(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
        fb(1,1) { ...
            mustBeNonsparse, mustBeNonempty, mustBeFloat, ...
            mustBeNonpositive ...
            };
    end
    
    methods
        
        function obj = Level(domain,exec,l)
            
            if ~nargin
                return;
            end
            
            narginchk(3,3)
            
            obj.Domain = tribosolver.internal.Domain(domain,exec,l);
            obj.Results.p = initP(obj.Domain.x,obj.Domain.y);
            obj.Results.x = obj.Domain.x;
            obj.Results.y = obj.Domain.y;
            obj.h = initH(obj.Domain.x,obj.Domain.y);
            
            [nx,ny] = size(obj.Domain.x);
            obj.k = initK(nx,ny,obj.Domain.dx,obj.Domain.dy);
            
            obj.fb = cast(-2*pi/3,"like",obj.h);
            obj.p_rhs = zeros(nx,ny,"like",obj.h);
            obj.p_old = zeros(nx,ny,"like",obj.h);
            
            padX = obj.Domain.nx*3-1;
            padY = obj.Domain.ny*3-1;
            obj.k_fft = fft2(obj.k,padX,padY);
            
            obj.calcDeformation();
            
        end
        
        function calcDeformation(obj)
            
            if numel(obj.k) < (32*32)
                
                obj.Results.w = conv2(obj.Results.p,obj.k,'same');
            else
                
                nx = obj.Domain.nx;
                ny = obj.Domain.ny;
                [padX,padY] = size(obj.k_fft);
                w=ifft2(fft2(obj.Results.p,padX,padY) .* obj.k_fft);
                obj.Results.w=w(nx:(2*nx-1),ny:(2*ny-1));
                
            end
            
            obj.Results.h = obj.Results.h0 + obj.h + obj.Results.w;
            
        end
        
    end
end

function p = initP(x,y)
value = x.^2 + y.^2;
value(value>1) = 1;
p=sqrt(1 - value);
end

function h = initH(x,y)
h = 0.5*x.^2 + 0.5*y.^2;
end

function K = initK(nx,ny,dx,dy)

x = cast( (0 : nx - 1)', "like", dx );
y = cast( (0 : ny - 1),  "like", dy );

x = (x + 0.5).*dx;
y = (y + 0.5).*dy;

xx = repmat(reshape(x,[],1),1,ny);
yy = repmat(reshape(y,1,[]),nx,1);

xm = xx - dx;
ym = yy - dy;

K = 2/(pi*pi)* ...
    (abs(xx).*asinh(yy ./ xx) ...
    + abs(yy).*asinh(xx ./ yy) ...
    - abs(xm).*asinh(yy ./ xm) ...
    - abs(yy).*asinh(xm ./ yy) ...
    - abs(xx).*asinh(ym ./ xx) ...
    - abs(ym).*asinh(xx ./ ym) ...
    + abs(xm).*asinh(ym ./ xm) ...
    + abs(ym).*asinh(xm ./ ym));

% Unboxing Kernel
DL = flip(K,2); DR = K; UL = flip(DL,1); UR = flip(K,1);
K=[UL(1:end-1,1:end-1) UR(1:end-1,1:end);DL(1:end,1:end-1) DR];

end