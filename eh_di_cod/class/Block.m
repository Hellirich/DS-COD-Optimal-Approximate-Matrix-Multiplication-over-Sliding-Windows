classdef Block
  properties
    mx
    my
    l
    t_start
    t_end
    X
    Y
    idx
    X_size
    Y_size
    sizex
    sizey
    size
  end

  methods
    function obj = Block(mx, my, l, t)
        obj.mx = mx;
        obj.my = my;
        obj.l = l;
        obj.t_start = t;
        obj.t_end = t;
        obj.X = zeros(mx, l);
        obj.Y = zeros(my, l);
        obj.idx = 0;
        obj.X_size = 0;
        obj.Y_size = 0;
        obj.sizex = 0;
        obj.sizey = 0;
        obj.size = 0;
    end

    function obj = clear(obj, t)
        obj.t_start = t;
        obj.t_end = t;
        obj.X = zeros(obj.mx, obj.l);
        obj.Y = zeros(obj.my, obj.l);
        obj.idx = 0;
        obj.sizex = 0;
        obj.sizey = 0;
        obj.size = 0;
    end

    function obj = insert(obj, x, y, t)
        obj.X(:, obj.idx + 1) = x;
        obj.Y(:, obj.idx + 1) = y;
        obj.sizex = obj.sizex + norm(x)^2;
        obj.sizey = obj.sizey + norm(y)^2;
        obj.size = sqrt(obj.sizex) * sqrt(obj.sizey);
        obj.idx = obj.idx + 1;
        obj.t_end = t;
    end
  end
end