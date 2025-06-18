classdef FastDumpSnapshotCod < handle
    
    properties
        N
        mx
        my
        l
        queueS % snapshots queue
        lenS % snapshot number
        sketchA % residual sketch
        sketchB
        K_A % ATA
        K_B % BTB
        idx
        threshold
        last_dump_time
        sketch_size % residul maximun sigular value
        level % threshold level for adaptive DS-COD
    end
    
    methods

        function obj = FastDumpSnapshotCod(N, mx, my, l, theta, t)
            obj.N = N;
            obj.mx = mx;
            obj.my = my;
            obj.l = l;
            obj.queueS = {};
            obj.lenS = 0;
            obj.sketchA = zeros(mx, 2*l);
            obj.sketchB = zeros(my, 2*l);
            obj.K_A = zeros(2*l, 2*l);
            obj.K_B = zeros(2*l, 2*l);
            obj.idx = 0;
            obj.threshold = theta;
            obj.last_dump_time = t;
            obj.sketch_size = 0;
            obj.level = 1;
        end

        function update(obj, x, y, current_time)

            norm_xy = norm(x)*norm(y);
            if norm_xy >= obj.threshold
                snap = snapshot(x,y,obj.last_dump_time+1,current_time);
                obj.queueS{end+1} = snap;
                obj.lenS = obj.lenS + snap.len;
                obj.last_dump_time = current_time;
                return;
            end

            obj.sketch_size = obj.sketch_size + norm_xy;
            if obj.idx + 1 >= 2 * obj.l % 2l buffer is full, cs to l
                obj.idx = obj.idx + 1;
                obj.sketchA(:, obj.idx) = x;
                obj.sketchB(:, obj.idx) = y;

                % time cost badly
                [obj.sketchA, obj.sketchB, obj.idx, sigma] = cs(obj.sketchA, obj.sketchB, obj.l);

                nonsnap_idx = find(sigma < obj.threshold, 1, 'first');
                if isempty(nonsnap_idx)
                    snapcount = length(sigma);
                else
                    snapcount = nonsnap_idx - 1;
                end
                % dump snapshot
                if snapcount > 0
                    a = obj.sketchA(:,1:snapcount);
                    b = obj.sketchB(:,1:snapcount);
                    snap = snapshot(a,b,obj.last_dump_time+1,current_time);
                    obj.queueS{end+1} = snap;
                    obj.lenS = obj.lenS + snap.len;
                    obj.sketchA(:,1:2*obj.l-snapcount) = obj.sketchA(:, snapcount+1:end);
                    obj.sketchB(:,1:2*obj.l-snapcount) = obj.sketchB(:, snapcount+1:end);
                    obj.sketchA(:, 2*obj.l-snapcount+1:end) = 0;
                    obj.sketchB(:, 2*obj.l-snapcount+1:end) = 0;
                    obj.idx = obj.idx - snapcount;
                    obj.last_dump_time = current_time;
                end

                % update K
                obj.K_A(:,:) = 0; obj.K_B(:,:) = 0;
                obj.K_A(1:obj.idx,1:obj.idx) = obj.sketchA(:,1:obj.idx)' * obj.sketchA(:,1:obj.idx);
                obj.K_B(1:obj.idx,1:obj.idx) = obj.sketchB(:,1:obj.idx)' * obj.sketchB(:,1:obj.idx);
                if snapcount < length(sigma)
                    obj.sketch_size = sigma(snapcount + 1);
                else
                    obj.sketch_size = 0;
                end

            else 
                % update sketch and K
                if obj.idx == 0
                    obj.K_A(1,1) = x' * x;
                    obj.K_B(1,1) = y' * y;
                else
                    obj.K_A(1:obj.idx, obj.idx+1) = sparse(obj.sketchA(:,1:obj.idx)') * x;
                    obj.K_A(obj.idx+1, 1:obj.idx) = obj.K_A(1:obj.idx, obj.idx+1)';
                    obj.K_A(obj.idx+1, obj.idx+1) = x' * x;
                    obj.K_B(1:obj.idx, obj.idx+1) = sparse(obj.sketchB(:,1:obj.idx)') * y;
                    obj.K_B(obj.idx+1, 1:obj.idx) = obj.K_B(1:obj.idx, obj.idx+1)';
                    obj.K_B(obj.idx+1, obj.idx+1) = y' * y;
                end
                obj.idx = obj.idx + 1;
                obj.sketchA(:, obj.idx) = x;
                obj.sketchB(:, obj.idx) = y;

                % dump snapshots
                if obj.sketch_size > obj.threshold
                    Ka = obj.K_A(1:obj.idx,1:obj.idx);
                    Kb = obj.K_B(1:obj.idx,1:obj.idx);
                    try
                        [L_a, D_a] = ldl(real(Ka));
                    catch
                        disp(Ka);
                    end
                    Rx=sqrtm(D_a)*L_a';
                    [L_b, D_b] = ldl(real(Kb));
                    Ry=sqrtm(D_b)*L_b';
                    [U, S, V] = svd(full(Rx * Ry'));

                    RxT = Rx';
                    RyT = Ry';
                    S = diag(S);
                    nonsnap_idx = find(S < obj.threshold, 1, 'first');
                    if isempty(nonsnap_idx)
                        snapcount = length(S);
                    else
                        snapcount = nonsnap_idx - 1;
                    end
                    for j = 1 : snapcount
                        a = sqrt(1/S(j)) * obj.sketchA(:,1:obj.idx) * RyT * V(:, j);
                        b = sqrt(1/S(j)) * obj.sketchB(:,1:obj.idx) * RxT * U(:, j);
                        snap = snapshot(a,b,obj.last_dump_time+1,current_time);
                        obj.queueS{end+1} = snap;
                        obj.lenS = obj.lenS + snap.len;
                        
                        % update sketch and K
                        obj.sketchA(:,1:obj.idx) = obj.sketchA(:,1:obj.idx) - 1/S(j) * obj.sketchA(:,1:obj.idx) * (RyT * V(:, j)) * (U(:, j)' * Rx);
                        obj.sketchB(:,1:obj.idx) = obj.sketchB(:,1:obj.idx) - 1/S(j) * obj.sketchB(:,1:obj.idx) * (RxT * U(:, j)) * (V(:, j)' * Ry);
                        ATA = Ka - (1/S(j)) * Ka * RyT * V(:, j) * U(:, j)' * Rx;
                        obj.K_A(1:obj.idx,1:obj.idx) = ATA - (1/S(j)) * RxT * U(:, j) * (V(:, j)' * Ry * ATA);
                        BTB = Kb - (1/S(j)) * Kb * RxT * U(:, j) * V(:, j)' * Ry;
                        obj.K_B(1:obj.idx,1:obj.idx) = BTB - (1/S(j)) * RyT * V(:, j) * (U(:, j)' * Rx * BTB);
                        % obj.K_A = obj.sketchA' * obj.sketchA;
                        % obj.K_B = obj.sketchB' * obj.sketchB;
                    end
                    if snapcount < length(S)
                        obj.sketch_size = S(snapcount + 1);
                    else
                        obj.sketchA(:,:) = 0;
                        obj.sketchB(:,:) = 0;
                        obj.idx = 0;
                        obj.K_A(:,:) = 0;
                        obj.K_B(:,:) = 0;
                        obj.sketch_size = 0;
                    end
                end
            end
        end

        function expire_S(obj, i)
            while ~isempty(obj.queueS) && obj.queueS{1}.t + obj.N < i
                obj.lenS = obj.lenS - obj.queueS{1}.len;
                obj.queueS(1) = [];
            end
        end

        function cap_S(obj)
            while ~isempty(obj.queueS) && obj.lenS > 3 * obj.l
                obj.lenS = obj.lenS - obj.queueS{1}.len;
                obj.queueS(1) = [];
            end
        end

        % adapt level for aDS
        function adapt_level(obj, L)
            obj.threshold = 2 ^ L * (obj.N / obj.l);
        end

        function [A,B] = query(obj)
            A = obj.sketchA(:,1:obj.idx);
            B = obj.sketchB(:,1:obj.idx);

            if ~isempty(obj.queueS)
                jdx = 1;
                Sx = zeros(obj.mx, 3*obj.l);
                Sy = zeros(obj.my, 3*obj.l);
                for j = 1:length(obj.queueS)
                    len = obj.queueS{j}.len;
                    Sx(:, jdx:jdx+len-1) = obj.queueS{j}.x;
                    Sy(:, jdx:jdx+len-1) = obj.queueS{j}.y; 
                    if jdx > 2*obj.l
                        [Sx,Sy,jdx] = cs(Sx,Sy,obj.l);
                    end
                    jdx = jdx + len;
                end
                A = [Sx, A];
                B = [Sy, B];
            end

        end
    
        function restart(obj, aux)
            obj.sketchA = aux.sketchA;
            obj.sketchB = aux.sketchB;
            obj.K_A = aux.K_A;
            obj.K_B = aux.K_B;
            obj.queueS = aux.queueS;
            obj.lenS = aux.lenS;
            obj.idx = aux.idx;
            obj.last_dump_time = aux.last_dump_time;
            obj.sketch_size = aux.sketch_size;
        end
    
    end
end

