classdef Domain < handle
    
    % Properties that define the global domain on which the problem is to
    % be solved. Internal solver use only.
    %
    % Copyright (C) 2021 theComputeKid
    
    properties(SetAccess=immutable)
        
        type(1,1) string { ...
            ismember(type,["constant","chebyshev"]) ...
            } = "constant";
        
        dx(:,1) { ...
            mustBeFloat, mustBeNonsparse, mustBeNonempty, ...
            mustBePositive ...
            } = 1;
        
        dy(:,1) { ...
            mustBeFloat, mustBeNonsparse, mustBeNonempty, ...
            mustBePositive ...
            } = 1;
        
        x(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
        y(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
    end
    
    methods
        
        function obj = Domain(domain,exec,l)
            
            if ~nargin
                return;
            end
            
            narginchk(3,3)
            
            p = exec.getProto();
            t = underlyingType(p);
            
            if strcmpi(domain.typeX,"constant")
                nx = cast(domain.nx,t);
                
                for i = 2:l
                    nx = nx*2 - 1;
                end
                
                obj.dx = (domain.xout - domain.xin)/(nx - 1);
                obj.dx = cast(obj.dx,t);
                obj.x = domain.xin + (0:nx-1)*obj.dx;
            end
            
            if strcmpi(domain.typeY,"constant")
                ny = cast(domain.ny,t);
                
                for i = 2:l
                    ny = ny*2 - 1;
                end
                
                obj.dy = (domain.yout - domain.yin)/(ny - 1);
                obj.dy = cast(obj.dy,t);
                obj.y = domain.yin + (0:ny-1)*obj.dy;
            end
            
            [obj.x,obj.y] = meshgrid(obj.x,obj.y);
            %             obj.dx = gradient(obj.x);
            %             obj.dy = gradient(obj.y);
            
        end
        
    end
end