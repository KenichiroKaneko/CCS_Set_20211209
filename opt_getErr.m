function [obj, err] = opt_getErr(obj, temp)

    % メモリーに同じ遺伝子が存在したら誤差を読み取る
    for i = 1:size(obj.gene_memory, 1)
        if temp.gene == obj.gene_memory(i,:)
            err = obj.err_memory(i);
            return
        end
    end

    % CCS
    % load('vars_verybig2033');
    adopt_BZ = [1, temp.gene(1:temp.BZ_gene_num), 1, fliplr(temp.gene(1:temp.BZ_gene_num))];
    adopt_FL = [1, temp.gene(temp.FL_gene_num+1:end), 1, fliplr(temp.gene(temp.FL_gene_num+1:end))];
    ind_BZ = find(adopt_BZ == 1);
    ind_FL = find(adopt_FL == 1);
    % refs = ["vars_verybig2033", "vars_verybig1000", "vars_verybig0800"];
    refs = ["vars_verybig2033", "UTST_numel_0720Ip50","UTST_numel_0760Ip50","UTST_numel_0800Ip50","UTST_numel_0840Ip50","UTST_numel_0880Ip50",...
            "UTST_numel_0920Ip50","UTST_numel_1000Ip50"];
    % refs = ["UTST_numel_2033", "UTST_numel_0720Ip50","UTST_numel_0760Ip50","UTST_numel_0800Ip50","UTST_numel_0840Ip50","UTST_numel_0880Ip50",...
    %         "UTST_numel_0920Ip50","UTST_numel_1000Ip50"];

    for i = 1:length(refs)

        load(['./vars_opt/'+ refs(i)], 'PARAM', 'CONFIG', 'FFout', 'SENSOR_TPRB', 'SENSOR_NPRB', 'SENSOR_FLXLP', 'CCSDAT', 'ExtCOIL', 'WALL');
        CONFIG.ShowFig = 0;
        % SENSORを作り直す
        SENSOR_TPRB.R = SENSOR_TPRB.R(ind_BZ);
        SENSOR_TPRB.Z = SENSOR_TPRB.Z(ind_BZ);
        SENSOR_TPRB.TET = SENSOR_TPRB.TET(ind_BZ);
        SENSOR_TPRB.TPRB = SENSOR_TPRB.TPRB(ind_BZ);
        SENSOR_TPRB.NUM = length(SENSOR_TPRB.TPRB);
        SENSOR_FLXLP.R = SENSOR_FLXLP.R(ind_FL);
        SENSOR_FLXLP.Z = SENSOR_FLXLP.Z(ind_FL);
        SENSOR_FLXLP.TET = SENSOR_FLXLP.TET(ind_FL);
        SENSOR_FLXLP.FLXLP = SENSOR_FLXLP.FLXLP(ind_FL);
        SENSOR_FLXLP.NUM = length(SENSOR_FLXLP.FLXLP);

        % 順問題の解行列FF、CCS点と各センサー及び自分以外のCCS点との関係式行列AAを作成
        FF = make_FF(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT);
        [FC, BR, BZ, PSIFLX, PSIC, AA, FF] = FORM(PARAM, CONFIG, FF, ...
        ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

        % 重み付け
        [FF, FC, AA, factors] = weight(FF,FC,AA,SENSOR_NPRB,SENSOR_TPRB,SENSOR_FLXLP,CCSDAT);

        % 逆問題を解く
        FFout_ = tsvd(PARAM, CONFIG, AA, FF, FC, ...
        SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);
        errs(i) = sum(abs(FFout - FFout_));
    end

    
    err = max(errs);

    % メモリーがobj.N*3個以下、errがほどほどに良かったら記憶する
    if size(obj.gene_memory, 1) < obj.N*3
        i = size(obj.gene_memory, 1) + 1;
        obj.gene_memory(i,:) = temp.gene;
        obj.err_memory(i) = err;
    elseif size(obj.gene_memory, 1) == obj.N*3
        max_err_me = obj.err_memory(1);
        for i = 2:size(obj.gene_memory, 1)
            if obj.err_memory(i) > max_err_me
                max_err_me = obj.err_memory(i);
            end
        end
        if err < 1.5 * max_err_me
            obj.gene_memory(1,:) = [];
            obj.err_memory(1) = [];
            obj.gene_memory(20,:) = temp.gene;
            obj.err_memory(20) = err;
        end
    else
        disp('something error happened at opt_getErr');
        i = size(obj.gene_memory, 1);
        obj.gene_memory(20:i,:) = [];
        obj.err_memory(20:i) = [];
    end


    % err = sum(abs(FF' - (AA * FFout)));

    % % 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
    % [psi, CCR, CCZ, DELGE, RCCS, ZCCS, FL2, BZ2, BR2] = INTER(PARAM, 0, FFout, ExtCOIL,...
    % SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % % 誤差評価する
    % err = evaluate_LCFS(psi, REF, PARAM, CCR, CCZ, 0);
end