classdef Results < handle
    
    % Holds result from a solution sweep
    %
    % Copyright (C) 2021 theComputeKid
    
    properties
        
        p(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
        h(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
        w(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
        h0(1,1) {mustBeFinite, mustBeNonsparse, mustBeFloat} = -0.3;
        
        x(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
        y(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
    end
    
    methods
        
        function plotH(obj)
            title("Film Thickness")
            surf(obj.x,obj.y,obj.h);
            xlabel("X"); ylabel("Y");zlabel("H");
        end
        
        function plotW(obj)
            title("Deformation")
            surf(obj.x,obj.y,obj.w);
            xlabel("X"); ylabel("Y");zlabel("H");
        end
        
        function plotP(obj)
            title("Fluid Pressure")
            surf(obj.x,obj.y,obj.p);
            xlabel("X"); ylabel("Y");zlabel("H");
        end
        
    end
end