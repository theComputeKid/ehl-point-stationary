classdef Results
    
    % Holds results from a solution sweep. Plot functions are provided.
    %
    % Copyright (C) 2021 theComputeKid
    
    properties
        
        % Fluid pressure solution
        p(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
        % Film thickness solution
        h(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
        % Deformation solution
        w(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
        % Final H0 constant value
        h0(1,1) {mustBeFinite, mustBeNonsparse, mustBeFloat} = -0.53;
        
        % Domain X values
        x(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
        % Domain Y values
        y(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
    end
    
    methods
        
        function plotH(obj)
            % Plot film thickness
            title("Film Thickness")
            surf(obj.x,obj.y,obj.h);
            xlabel("X"); ylabel("Y");zlabel("H");
        end
        
        function plotW(obj)
            % Plot deformation solution
            title("Deformation")
            surf(obj.x,obj.y,obj.w);
            xlabel("X"); ylabel("Y");zlabel("W");
        end
        
        function plotP(obj)
            % PLot fluid pressure solution
            title("Fluid Pressure")
            surf(obj.x,obj.y,obj.p);
            xlabel("X"); ylabel("Y");zlabel("P");
        end
        
    end
end