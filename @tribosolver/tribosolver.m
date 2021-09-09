classdef tribosolver < handle
    % Container for a Point Contact ElastoHydroDynamic Lubrication Model.
    % Parameters can be provided in any order. The Domain and Moes
    % parameters are required, the other arguments are optional.
    %
    % Syntax:
    % model = tribosolver(Domain,Moes,Execution,Relaxation,Tolerance)
    %
    % Copyright (C) 2021 theComputeKid
    
    properties(SetAccess=immutable)
        
        Domain(1,1) tribosolver.Domain
        
        Moes(1,1) tribosolver.Moes
        
        Tolerance(1,1) tribosolver.Tolerance
        
        Relaxation(1,1) tribosolver.Relaxation
        
        Execution(1,1) tribosolver.Execution
        
    end
    
    properties
        
        Levels(:,1) tribosolver.internal.Level
        
    end
    
    methods
        
        function obj = tribosolver(varargin)
            
            % Constructor for the model
            %
            % Syntax:
            % obj = tribosolver(domain,moes)
            % obj = tribosolver(domain,moes,exec)
            % obj = tribosolver(domain,moes,exec,relax)
            % obj = tribosolver(domain,moes,exec,relax,tol)
            %
            % Copyright (C) 2021 theComputeKid
            
            % Moes and Domain Parameters are necessary, the others are
            % optional.
            for i = 1:nargin
                
                if isa(varargin{i},"tribosolver.Domain")
                    obj.Domain = varargin{i};
                    hasSetDomain = true;
                elseif isa(varargin{i},"tribosolver.Moes")
                    obj.Moes = varargin{i};
                    hasSetMoes = true;
                elseif isa(varargin{i},"tribosolver.Tolerance")
                    obj.Tolerance = varargin{i};
                elseif isa(varargin{i},"tribosolver.Relaxation")
                    obj.Relaxation = varargin{i};
                elseif isa(varargin{i},"tribosolver.Execution")
                    obj.Execution = varargin{i};
                end
                
            end
            
            assert(hasSetDomain, ...
                "Tribosolver: Domain parameters not provided")
            
            assert(hasSetMoes, ...
                "Tribosolver: Moes parameters not provided")
            
        end
        
        results = solve(obj,verbosity);
        
    end
    
    methods(Access=private)
        
        initLevels(obj);
        mgFull(obj,l);
        mgCycle(obj,l);
        relaxP(obj,l);
        relaxH(obj,l);
        coarsen(obj,l);
        coarsenP(obj,l);
        coarsenPRHS(obj,l);
        coarsenFB(obj,l);
        refine(obj,l);
        mgFullInterp(obj,l);
        
    end
end