classdef Moes
    
    % Non-dimensional Moes Parameters used in the model.
    %
    % Syntax:
    % obj = Execution(basePrecision,device, verbosity)
    %
    % Input arguments:
    % M: Moes M parameter
    % L: Moes L parameter
    % H0: Initial H0 solution guess
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
        
        H0(1,1) { ...
            mustBeFinite, mustBeNonsparse, mustBeFloat ...
            } = -0.3;
    end
    
    methods
        
        function obj = Moes(M,L,H0)
            
            % Syntax:
            % obj = Moes(M)
            % obj = Moes(M,L)
            
            if ~nargin
                return;
            end
            
            narginchk(1,3)
            
            obj.M = M;
            
            if nargin > 1
                obj.L = L;
            end
            
            if nargin > 2
                obj.H0 = H0;
            end
        end
        
    end
    
end
