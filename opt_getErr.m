function err = opt_getErr(obj, temp)
    % CCS
    load('vars_verybig2033');
    adopt_BZ = [1, temp.gene(1:temp.BZ_gene_num), 1, fliplr(temp.gene(1:temp.BZ_gene_num))];
    adopt_FL = [1, temp.gene(temp.BZ_gene_num+1:end), 1, fliplr(temp.gene(temp.BZ_gene_num+1:end))];
    ind_BZ = find(adopt_BZ == 1);
    ind_FL = find(adopt_FL == 1);
    refs = ["vars_verybig2033", "vars_verybig1000", "vars_verybig0800"];

    for i = 1:3

        load(refs(i), 'PARAM', 'CONFIG', 'FFout', 'SENSOR_TPRB', 'SENSOR_NPRB', 'SENSOR_FLXLP', 'CCSDAT', 'ExtCOIL', 'WALL');
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
    % err = sum(abs(FF' - (AA * FFout)));

    % % 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
    % [psi, CCR, CCZ, DELGE, RCCS, ZCCS, FL2, BZ2, BR2] = INTER(PARAM, 0, FFout, ExtCOIL,...
    % SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % % 誤差評価する
    % err = evaluate_LCFS(psi, REF, PARAM, CCR, CCZ, 0);
end