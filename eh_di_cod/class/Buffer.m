classdef Buffer

    properties
        X_buffer
        Y_buffer
        X_size
        Y_size
        X_nnz
        Y_nnz
        X_entry_start
        Y_entry_start
        columns
        mx
        my
        t_start
        t_end
    end
    
    methods
        function obj = Buffer(mx, my)
            obj.X_size = 0;
            obj.Y_size = 0;
            obj.X_nnz = 0;
            obj.Y_nnz = 0;
            obj.X_entry_start = 1;
            obj.Y_entry_start = 1;
            obj.columns = 0;
            obj.mx = mx;
            obj.my = my;
            obj.t_start = 1;
            obj.t_end = 1;
        end
        
        function obj = reset(obj, t)
            obj.X_size = 0;
            obj.Y_size = 0;
            obj.X_nnz = 0;
            obj.Y_nnz = 0;
            obj.X_entry_start = obj.X_entry_start + obj.X_nnz;
            obj.Y_entry_start = obj.Y_entry_start + obj.Y_nnz;
            obj.columns = 0;
            obj.t_start = t;
            obj.t_end = t;
        end
    end
end

