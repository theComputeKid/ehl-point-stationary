classdef Results < handle
    
    % Holds result from a solution sweep
    %
    % Copyright (C) 2021 theComputeKid
    
    properties
        
        p(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
        w(:,:) {mustBeNonsparse, mustBeNonempty, mustBeFloat} = 0;
        
        h0(1,1) {mustBeFinite, mustBeNonsparse, mustBeFloat} = -0.3;
        
    end
    
end