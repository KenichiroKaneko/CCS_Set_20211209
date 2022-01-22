classdef opt_GA < handle
    properties (Constant)
        % 個体数
        N = 20;
        % 世代数
        itr = 4000;
        % 許容誤差
        EPS = 0.001;
    end

    properties
        pool = opt_GA.init(opt_GA.N);
        pool_next = opt_GA.init(opt_GA.N);
        err = 1;
        min_err = 2^16;
        goodgene = [];
        gooditr = 1;
        gene_memory = [];
        err_memory = [];
    end

    methods (Static)
        function pool = init(N)
            for i = 1:N
                % Matlabはゴミ
                pool(i) = opt_Chromosome();

                % x = zeros(1,pool(1).BZ_gene_num);
                % y = zeros(1,pool(1).FL_gene_num);
                % ind_x = [1  2  4  5  7 11 13 16 19 22 25 27 31 41 42 45 46 50 52];
                % ind_y = [3  4  6 11 17 20 27 30 31 33 36 37 42 44 45 46 47 48 52];  
                % x(ind_x) = 1;
                % y(ind_y) = 1;

                x = pool(i).gene(1:pool(i).BZ_gene_num);
                y = pool(i).gene(pool(i).BZ_gene_num+1:end);
                x = x(randperm(length(x)));
                y = y(randperm(length(y)));
                pool(i).gene = [x y];
                obj.gene_memory(i,:) = pool(i).gene;
            end
        end
    end

    methods
        function evolve(obj)
            for i = 1:obj.itr
                obj.getVal();
                obj.printStatus(i);
                obj.tournamentSelection();
                % obj.rouletteSelection();
                obj.crossOverAll();
                obj.crossChangeAll();
                obj.saveStatus(i);
                obj.mutateAll();
                if i > 10
                    if sum(obj.err(end-1:end)) < obj.EPS
                        disp(['GA finished at itr:' num2str(i)]);
                        break
                    end
                end
                for j = 1:obj.N
                    if obj.pool(j).err < obj.min_err
                        obj.min_err = obj.pool(j).err;
                        obj.goodgene = obj.pool(j).gene;
                        obj.gooditr = i;
                    end
                end
            end
        end

        function crossOverAll(obj)
            % 交叉
            for i = 1:obj.N
                rand1 = randi(obj.N);
                [obj.pool(i), obj.pool(rand1)] = obj.pool(i).crossOver(i, rand1, obj.pool(rand1));
            end
        end

        function crossChangeAll(obj)
            for i = 1:obj.N
                obj.pool(i) = obj.pool(i).crossChange();
            end
        end

        function obj = mutateAll(obj)
            % 突然変異
            for i = 1:obj.N
                obj.pool(i) = obj.pool(i).mutate();
            end
        end

        function obj = rouletteSelection(obj)
            % ルーレット選択
            C = 0;
            errs = [];
            for i = 1:obj.N
                errs(i) = obj.pool(i).err;
                C = C + errs(i);
            end
            errs = (C - errs).^2;
            C = sum(errs);
            rates = (errs)./C;
            roulette = rates(1);
            for i = 1:obj.N-1
                roulette(i+1) = roulette(i) + rates(i+1);
            end
            i = 1;
            while i <= obj.N
                r = rand;
                for j = 1:obj.N
                    if r <= roulette(j)
                        obj.pool_next(i) = copy(obj.pool(j));
                        break
                    end
                end
                i = i+1;
            end
            obj.pool = obj.pool_next;
        end

        function obj = tournamentSelection(obj)
            % トーナメント選択
            i = 1;
            while i <= obj.N
                box = 2;
                win = 1;
                rands = [];
                offspring_errs = [];
                for j = 1:box
                    rands(j) = randi(obj.N);
                    offspring_errs(j) = obj.pool(rands(j)).err;
                end
                for j = 1:win
                    [M I] = min(offspring_errs);
                    offspring = copy(obj.pool(rands(I)));
                    offspring.parentX = rands(I);
                    offspring.parentY = 0;
                    offspring.M = ' ';
                    offspring.X = ' ';
                    obj.pool_next(i) = offspring;
                    i = i + 1;
                    offspring_errs(I) = [];
                    rands(I) = [];
                end
            end
            obj.pool = obj.pool_next;
        end

        function obj = getVal(obj)
            for i = 1:obj.N
                temp = copy(obj.pool(i));
                [obj, obj.pool(i).err] = opt_getErr(obj, temp);
            end
        end

        function [obj, err] = opt_getErr2(obj, temp)
            BZ_i = find(temp.gene(1 : temp.BZ_gene_num) == 1);
            FL_i = find(temp.gene(temp.BZ_gene_num+1 : end) == 1);
            err1 = (abs(333 - sum([BZ_i FL_i])));
            if length([BZ_i FL_i]) < 10
                err2 = 100;
            else
                err2 = length([BZ_i FL_i]) - 10;
            end
            err = err1 + err2;
        end

        function printStatus(obj, i)
            if (rem(i, 1) == 0 || i == 1)
                % Log = fopen("optLog.txt", "a");
                % fprintf(Log, sprintf('generation : %03d\n', i));
                disp(sprintf('generation : %03d', i));
                disp(sprintf('gene / err / BZpos / FLpos'));
                errs = [];
                for j = 1:obj.N
                    errs(j) = obj.pool(j).err;
                end
                [errs I] = sort(errs);
                for j = 1:obj.N
                    k = I(j);
                    parents = sprintf('%02d:%02d-%02d%s',[k, obj.pool(k).parentX, obj.pool(k).parentY, obj.pool(k).M, obj.pool(k).X]);
                    BZ_i = find(obj.pool(k).gene(1 : obj.pool(k).BZ_gene_num) == 1);
                    BZ_t = sprintf(repmat('% 3d', [1, length(BZ_i)]), BZ_i);
                    FL_i = find(obj.pool(k).gene(obj.pool(k).BZ_gene_num+1 : end) == 1);
                    FL_t = sprintf(repmat('% 3d', [1, length(FL_i)]), FL_i);
                    err_t = sprintf(' %.5e', obj.pool(k).err);
                    disp([parents err_t BZ_t FL_t]);
                    % fprintf(Log, [parents err_t BZ_t FL_t '\n']);
                end
            end
        end

        function saveStatus(obj, i)
            % if rem(i, 10) == 0
            %     save(['vars_opt_generation_', num2str(i)]);
            % end
            min_err = obj.pool(1).err;
            for j = 2:obj.N
                min_err = min(min_err, obj.pool(j).err);
            end
            obj.err(i) = min_err;
        end
    end
end