classdef Execution
    
    % Execution strategy of the solver.
    %
    % Syntax:
    % obj = Execution(basePrecision,device, verbosity)
    %
    % Input arguments:
    % BasePrec: Precision of the calculations
    % Device: Choose whether to run on the CPU or GPU. Caution: GPU is not
    % optimized yet.
    % Verbosity: 1 for text update, 2 for graphical plot update during the
    % solution process (slows down solution process)
    %
    % Copyright (C) 2021 theComputeKid
    
    properties(SetAccess=immutable)
        
        BasePrecision(1,1) string { ...
            ismember(BasePrecision,["double","single"]) ...
            } = "double";
        
        Device(1,1) string { ...
            ismember(Device,["gpu","cpu_seq","cpu_par","gpu"]) ...
            } = "cpu_seq";
        
        Verbosity(1,1) uint64 {mustBeNonempty} = false;
        
    end
    
    methods
        
        function obj = Execution(basePrecision,device,verbosity)
            
            % Syntax:
            % obj = Execution(basePrecision,device)
            
            if ~nargin
                return;
            end
            
            narginchk(2,3)
            
            obj.BasePrecision = basePrecision;
            obj.Device = device;
            
            if (nargin > 2)
                obj.Verbosity = verbosity;
            end
            
        end
        
        function proto = getProto(obj)
            proto = cast([],obj.BasePrecision);
            if strcmpi(obj.Device,"gpu")
                proto = gpuArray(proto);
            end
        end
        
    end
end