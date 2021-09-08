classdef Execution
    
    % Execution strategy of the solver.
    %
    % Copyright (C) 2021 theComputeKid
    
    properties(SetAccess=immutable)
        
        BasePrecision(1,1) string { ...
            ismember(BasePrecision,["double","single"]) ...
            } = "double";
        
        Device(1,1) string { ...
            ismember(Device,["gpu","cpu"]) ...
            } = "cpu";
        
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