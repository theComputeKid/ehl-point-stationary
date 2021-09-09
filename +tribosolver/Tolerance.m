classdef Tolerance
    
    % Properties that define the tolerances for determining convergence of
    % the solution. Convergence is said to be achieved if the difference of
    % the old and new solution is within the tolerance specified here.
    %
    % Syntax:
    % obj = Tolerance(p,h,fb)
    %
    % Input arguments:
    % p: Fluid pressure
    % h: Film Thickness
    % fb: Force Balance
    %
    % Copyright (C) 2021 theComputeKid
    
    properties(SetAccess=immutable)
        
        p(1,1) double { ...
            mustBeFinite, mustBeNonsparse, mustBePositive ...
            } = 1e-5;
        
        h(1,1) double { ...
            mustBeFinite, mustBeNonsparse, mustBePositive ...
            } = 1e-5;
        
        fb(1,1) double { ...
            mustBeFinite, mustBeNonsparse, mustBePositive ...
            } = 1e-5;
        
    end
    
    methods
        
        function obj = Tolerance(p,h,fb)
            
            % Syntax:
            % obj = Tolerance(p,h,fb)
            
            if ~nargin
                return;
            end
            
            narginchk(3,3)
            
            obj.p = p;
            obj.h = h;
            obj.fb = fb;
            
        end
    end
end