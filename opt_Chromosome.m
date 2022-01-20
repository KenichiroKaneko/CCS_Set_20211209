% 最適化のループの中
classdef opt_Chromosome < handle & matlab.mixin.Copyable
    properties (Constant = true)
        % センサーの数
        BZ_num = 40;
        FL_num = 40;
        % z=0の点を固定＆上下対照にするときの遺伝子が1の数
        BZ_1gene_num = 18;
        FL_1gene_num = 18;
        % センサー候補地の数
        BZ_points = 106;
        FL_points = 106;
        % 上半分のセンサー候補地の数（遺伝子が0か1の数）
        BZ_gene_num = 52;
        FL_gene_num = 52;
        % % センサーの数
        % BZ_num = 5;
        % FL_num = 5;
        % % z=0の点を固定＆上下対照にするときの遺伝子が1の数
        % BZ_1gene_num = 5;
        % FL_1gene_num = 5;
        % % センサー候補地の数
        % BZ_points = 10;
        % FL_points = 10;
        % % 上半分のセンサー候補地の数（遺伝子が0か1の数）
        % BZ_gene_num = 10;
        % FL_gene_num = 10;
        %  遺伝子数
        L = opt_Chromosome.BZ_gene_num + opt_Chromosome.FL_gene_num;
        % 突然変異率
        M_rate = 0.005;
        % 交叉率
        C_rate = 0.1;
        % 特殊交叉率
        X_rate = 0.001;
    end

    properties (Constant = false)
        gene = opt_Chromosome.make_random_instance();
        err = 100;
        parentX = 0;
        parentY = 0;
        M = ' ';
        X = ' ';
    end

    methods (Static)
        function gene = make_random_instance()
            % ランダムな特徴を持ったインスタンスを生成
            x = [ones(1, opt_Chromosome.BZ_1gene_num), ... 
                zeros(1, opt_Chromosome.BZ_gene_num - opt_Chromosome.BZ_1gene_num)];
            y = [ones(1, opt_Chromosome.FL_1gene_num), ... 
                zeros(1, opt_Chromosome.FL_gene_num - opt_Chromosome.FL_1gene_num)];
            x_rand = x(randperm(length(x)));
            y_rand = y(randperm(numel(y)));
            x = zeros(1,opt_Chromosome.BZ_gene_num);
            y = zeros(1,opt_Chromosome.FL_gene_num);
            ind_x = [2 3:3:51 52];
            ind_y = [2 3:3:51 52];
            x(ind_x) = 1;
            y(ind_y) = 1;
            gene = [x, y];
        end
    end

    methods 
        function obj = mutate(obj)
            % 突然変異
            for i = 1:obj.L
                if rand <= obj.M_rate
                    obj.gene(i) = 1 - obj.gene(i);
                    obj.M = '*';
                end
            end
            % 突然変異によってセンサー数が増減したらランダムで元に戻す
            % 多い染色体を減らす
            % 1 染色体が40以上にならないようにする
            range = 1 : obj.BZ_gene_num;
            if sum(obj.gene(range)) ~= obj.BZ_1gene_num
                x = obj.BZ_1gene_num - sum(obj.gene(range));
                if x < 0
                    gene_type = 1;
                elseif x > 0
                    gene_type = 0;
                end
                obj = chromo_mdfy(obj, gene_type, abs(x), range, 'BZ');
            end
            range = obj.BZ_gene_num+1 : obj.L;
            if sum(obj.gene(range)) ~= obj.FL_1gene_num
                x = obj.FL_1gene_num - sum(obj.gene(range));
                if x < 0
                    gene_type = 1;
                elseif x > 0
                    gene_type = 0;
                end
                obj = chromo_mdfy(obj, gene_type, abs(x), range, 'FL');
            end
        end

        function [obj other] = chromo_mdfy(obj, gene_type, num, range, sen_type)
            % 突然変異の調整
            % gene_type の染色体を減らす
            i = find(obj.gene(range) == gene_type);
            i = i(randperm(numel(i)));
            i = i(1:num);
            if sen_type == 'FL'
                i = i + obj.BZ_gene_num;
            end
            obj.gene(i) = 1 - gene_type;
        end

        function [obj, other] = crossOver(obj, X, Y, other)
            % 交叉
            if 0 %rand <= obj.C_rate
                if rand < 0.5
                    % BZセンサーをYのものにする
                    range = 1:obj.BZ_gene_num;
                    temp = obj.gene(range);
                    obj.gene(range) = other.gene(range);
                    other.gene(range) = temp;
                    obj.parentX = X; obj.parentY = Y;
                    other.parentX = Y; other.parentY = X;
                else
                    % FLセンサーをYのものにする
                    range = obj.BZ_gene_num+1:obj.L;
                    temp = obj.gene(range);
                    obj.gene(range) = other.gene(range);
                    other.gene(range) = temp;
                    obj.parentX = Y; obj.parentY = X;
                    other.parentX = X; other.parentY = Y;
                end
            end
            if rand <= obj.C_rate
                i = find(rand(1, obj.L) < 0.5);
                temp = obj.gene;
                obj.gene(i) = other.gene(i);
                other.gene(i) = temp(i);
                obj.parentX = X; obj.parentY = Y;
                other.parentX = X; other.parentY = Y;
            end
        end

        function obj = crossChange(obj)
            % センサーの配置をBZとFLで入れ替える交叉
            if rand <= obj.X_rate
                tempBZ = obj.gene(1:obj.BZ_gene_num);
                tempFL = obj.gene(obj.BZ_gene_num+1:obj.L);
                obj.gene = [tempFL tempBZ];
                obj.X = 'X';
            end
        end
    end
end
