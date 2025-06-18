classdef snapshot
    properties
        x
        y
        len
        s
        t
    end
    
    methods
        function obj = snapshot(x,y,s,t)
            obj.x = x;
            obj.y = y;
            obj.len = size(x,2);
            obj.s = s;
            obj.t = t;
        end
    end
end

