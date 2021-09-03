classdef Domain
    
    % Properties that define the global domain on which the problem is to
    % be solved.
    %
    % Copyright (C) 2021 theComputeKid
    
    properties(SetAccess=immutable)
        
        mgl(1,1) uint64 { ...
            mustBeFinite, mustBeNonsparse, mustBeNonempty, ...
            mustBePositive ...
            } = 3;
        
        nx(1,1) uint64 { ...
            mustBeFinite, mustBeNonsparse, mustBeNonempty, ...
            mustBeGreaterThanOrEqual(nx,32) ...
            } = 32;
        
        ny(1,1) uint64 { ...
            mustBeFinite, mustBeNonsparse, mustBeNonempty, ...
            mustBeGreaterThanOrEqual(ny,32) ...
            } = 32;
        
        xin(1,1) double { ...
            mustBeFinite, mustBeNonsparse, mustBeNonempty, ...
            mustBeNegative ...
            } = -3;
        
        xout(1,1) double { ...
            mustBeFinite, mustBeNonsparse, mustBeNonempty, ...
            mustBePositive ...
            } = 3;
        
        yin(1,1) double { ...
            mustBeFinite, mustBeNonsparse, mustBeNonempty, ...
            mustBeNegative ...
            } = -3;
        
        yout(1,1) double { ...
            mustBeFinite, mustBeNonsparse, mustBeNonempty, ...
            mustBePositive ...
            } = 3;
        
        typeX(1,1) string { ...
            ismember(typeX,["constant","chebyshev"]) ...
            } = "constant";
        
        typeY(1,1) string { ...
            ismember(typeY,["constant","chebyshev"]) ...
            } = "constant";
    end
    
    methods
        
        function obj = Domain(xin,xout,nx,yin,yout,ny,mgl,typeX,typeY)
            
            % Syntax:
            % obj = Domain(xin,xout,nx,yin,yout,ny,mgl)
            % obj = Domain(xin,xout,nx,yin,yout,ny,mgl,typeX)
            % obj = Domain(xin,xout,nx,yin,yout,ny,mgl,typeX,typeY)
            
            if ~nargin
                return;
            end
            
            narginchk(7,9)
            
            obj.xin = xin;
            obj.xout = xout;
            obj.nx = nx;
            
            obj.yin = yin;
            obj.yout = yout;
            obj.ny = ny;
            
            obj.mgl = mgl;
            
            if nargin > 7
                obj.typeX = typeX;
            end
            
            if nargin > 8
                obj.typeY = typeY;
            end
            
        end
        
    end
end