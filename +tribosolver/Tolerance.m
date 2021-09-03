classdef Tolerance
    
    % Properties that define the tolerances for determining convergence of
    % the solution.
    %
    % Copyright (C) 2021 theComputeKid
    
    properties(SetAccess=immutable)
        
        p(1,1) double { ...
            mustBeFinite, mustBeNonsparse, mustBePositive ...
            } = 1e-6;
        
        h(1,1) double { ...
            mustBeFinite, mustBeNonsparse, mustBePositive ...
            } = 1e-6;
        
        fb(1,1) double { ...
            mustBeFinite, mustBeNonsparse, mustBePositive ...
            } = 1e-6;
        
    end
    
    methods
        
        function obj = Tolerance(p,h,fb)
            
            % Syntax:
            % obj = Tolerance(p,h)
            
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