classdef tribosolver < handle
    % TRIBOSOLVER:
    % Container for a Point Contact ElastoHydroDynamic Lubrication Model
    %
    % Syntax:
    % model = tribosolver(domain,moes)
    % model = tribosolver(domain,moes,exec)
    % model = tribosolver(domain,moes,exec,relax)
    % model = tribosolver(domain,moes,exec,relax,tol)
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
        
        function obj = tribosolver(domain,moes,exec,relax,tol)
            
            % Constructor for the model
            %
            % Syntax:
            % obj = tribosolver(domain,moes)
            % obj = tribosolver(domain,moes,exec)
            % obj = tribosolver(domain,moes,exec,relax)
            % obj = tribosolver(domain,moes,exec,relax,tol)
            %
            % Copyright (C) 2021 theComputeKid
            
            narginchk(2,5)
            
            obj.Moes = moes;
            obj.Domain = domain;
            
            if nargin > 2
                obj.Execution = exec;
            end
            
            if nargin > 3
                obj.Relaxation = relax;
            end
            
            if nargin > 4
                obj.Tolerance = tol;
            end
            
        end
        
        results = solve(obj,verbosity);
        
    end
    
    methods(Access=private)
        
        initLevels(obj);
        mgFull(obj,l);
        mgCycle(obj,l);
        relax(obj,l);
        coarsen(obj,l);
        refine(obj,l);
        mgFullInterp(obj,l);
        
    end
end