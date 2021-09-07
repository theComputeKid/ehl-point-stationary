classdef Relaxation

    % Relaxation properties used in the model.
    %
    % Copyright (C) 2021 theComputeKid

    properties

        switchEpsValue(1,1) double { ...
            mustBeFinite, mustBeNonsparse, mustBePositive ...
            } = 0.3;

        jac(1,1) double { ...
            mustBeFinite, mustBeNonsparse, mustBePositive ...
            } = 0.2;

        gs(1,1) double { ...
            mustBeFinite, mustBeNonsparse, mustBePositive ...
            } = 0.3;

        h0(1,1) double { ...
            mustBeFinite, mustBeNonsparse, mustBePositive ...
            } = 0.01;

        numCycles(1,1) { ...
            mustBeFinite, mustBeNonsparse, mustBeNonempty, ...
            mustBePositive ...
            } = 50;

        itPre(1,1) uint64 { ...
            mustBeFinite, mustBeNonsparse, mustBeNonempty, ...
            mustBePositive ...
            } = 3;

        itMain(1,1) uint64 { ...
            mustBeFinite, mustBeNonsparse, mustBeNonempty, ...
            mustBePositive ...
            } = 10;

        itPost(1,1) uint64 { ...
            mustBeFinite, mustBeNonsparse, mustBeNonempty, ...
            mustBePositive ...
            } = 3;

        gamma(1,1) uint64 { ...
            mustBeFinite, mustBeNonsparse, mustBeNonempty, ...
            mustBePositive ...
            } = 2;
    end

    methods

        function obj = Relaxation(jac,gs,h0,numCycles,gamma,switchEpsValue)

            % Syntax:
            % obj = Relaxation(jac,gs,h0)
            % obj = Relaxation(jac,gs,h0,switchEpsValue)

            if ~nargin
                return;
            end

            narginchk(3,6)

            obj.jac = jac;
            obj.gs = gs;
            obj.h0 = h0;

            if nargin > 3
                obj.numCycles = numCycles;
            end

            if nargin > 4
                obj.gamma = gamma;
            end

            if nargin > 5
                obj.switchEpsValue = switchEpsValue;
            end

        end

    end

end
