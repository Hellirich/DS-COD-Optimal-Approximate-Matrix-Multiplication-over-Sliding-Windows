classdef Candidate
    properties
        x
        y
        t
        weight
        rank = 1; % Default value
    end
    
    methods
        function obj = Candidate(x, y, t, weight)
            obj.x = x;
            obj.y = y;
            obj.t = t;
            obj.weight = weight;
            % rank is already initialized to 1 by default
        end
    end
end
