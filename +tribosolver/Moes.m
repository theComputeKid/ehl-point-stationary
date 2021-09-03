classdef Moes
    
    % Non-dimensional Moes Parameters used in the model.
    %
    % Copyright (C) 2021 theComputeKid
    
    properties(SetAccess=immutable)
        
        M(1,1) { ...
            mustBeFinite, mustBeNonsparse, mustBePositive, ...
            mustBeFloat ...
            } = 10;
        
        L(1,1) { ...
            mustBeFinite, mustBeNonsparse, mustBeNonnegative, ...
            mustBeFloat ...
            } = 0;
        
    end
    
    methods
        
        function obj = Moes(M,L)
            
            % Syntax:
            % obj = Moes(M)
            % obj = Moes(M,L)
            
            if ~nargin
                return;
            end
            
            narginchk(1,2)
            
            obj.M = M;
            
            if nargin > 1
                obj.L = L;
            end
            
        end
        
    end
    
end
