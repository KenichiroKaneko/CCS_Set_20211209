function CCS_UTST()
    %
    close all; clear all;
    type = 'exp';
    type = 'sol';
    [PARAM, CONFIG] = define_params(type);

    if CONFIG.Reload
        % データを最初から読み込む
        if CONFIG.DataType == 'exp'
            [FL, BZ, ExtCOIL, Ip] = load_lizzie(PARAM, CONFIG);
            POS = 0; REF = 0;
            save([PARAM.output_file_directory '/' PARAM.shotnum PARAM.date],...
            'CONFIG', 'FL', 'BZ', 'ExtCOIL', 'Ip', 'POS', 'REF');
            t = round((str2double(PARAM.time_CCS)) / 0.5);
            FL = FL(:, t)'; BZ = BZ(:, t)'; BR = zeros(1, length(BZ));
            ExtCOIL.I(3:10) = ExtCOIL.I_sig(1:8, t);

            % % 2Dプローブからreferenceのデータを作成
            % REF = load_2D_probe(PARAM);
            % hahaha;

        elseif CONFIG.DataType == 'sol'
            [FL, BZ, BR, ExtCOIL, REF, JEDDY, POS] = load_num_sol(PARAM, CONFIG); % OK
            t = 0; Ip = 0;
            save([PARAM.output_file_directory PARAM.num_sol_name],...
            'CONFIG', 'FL', 'BZ', 'BR', 'ExtCOIL', 'REF', 'JEDDY', 'POS');
        end
    else
        % 一度読み込んだデータをmatファイルから読み込む（はやい）
        % 入力はそのままで、PARAMだけいじりたい時はこちらを使う
        if CONFIG.DataType == 'exp'
            load([PARAM.output_file_directory '/' PARAM.shotnum PARAM.date]);
            t = round((str2double(PARAM.time_CCS)) / 0.5);
            FL = FL(:, t)'; BZ = BZ(:, t)'; BR = zeros(1, length(BZ));
            ExtCOIL.I(3:10) = ExtCOIL.I_sig(1:8, t);
        elseif CONFIG.DataType == 'sol'
            t = 0; Ip = 0;
            load([PARAM.output_file_directory PARAM.num_sol_name]);
        end
    end

    % 各パラメータを設定
    [PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP] = make_sensor(PARAM, CONFIG, BZ, BR, FL, POS);
    CCSDAT = make_CCS(PARAM);
    FF = make_FF(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT);
    % WALL = load_wall();
    WALL = loadwalldata2(PARAM, CONFIG);

    % 順問題の解行列FF、CCS点と各センサー及び自分以外のCCS点との関係式行列AAを作成
    [FC, BR, BZ, PSIFLX, PSIC, AA, FF] = FORM(PARAM, CONFIG, FF, ...
        ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % 重み付け
    [FF, FC, AA, factors] = weight(FF,FC,AA,SENSOR_NPRB,SENSOR_TPRB,SENSOR_FLXLP,CCSDAT);

    % 逆問題を解く
    FFout = tsvd(PARAM, CONFIG, AA, FF, FC, ...
        SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % プロット用に渦電流を計算する
    DISF = EDDYP(FFout, PARAM, CONFIG, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
    [psi, CCR, CCZ, DELGE, RCCS, ZCCS, FL2, BZ2, BR2] = INTER(PARAM, 0, FFout, ExtCOIL,...
    SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % 結果を表示
    if CONFIG.ShowFig
        show_results(PARAM, CONFIG, CCR, CCZ, CCSDAT, REF, psi,...
            SENSOR_TPRB, SENSOR_FLXLP, WALL,...
            ExtCOIL, DISF, t, Ip, FFout)
    end

    err_LCFS = evaluate_LCFS(psi, REF, PARAM, CCR, CCZ, 0)
    error('error description', A1)
    save('vars_sol_2033')
    % 
    % 解からもう一度逆問題を解く
    FF = AA * FFout; FF = FF';
    PARAM.dead_FL = []; PARAM.dead_BZ = [];
    PARAM.KUP = min(size(AA));
    
    % 重み付け
    [FF, FC, AA, factors] = weight(FF,FC,AA,SENSOR_NPRB,SENSOR_TPRB,SENSOR_FLXLP,CCSDAT);

    % 逆問題を解く
    FFout = tsvd(PARAM, CONFIG, AA, FF, FC, ...
        SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);
    
    % プロット用に渦電流を計算する
    DISF = EDDYP(FFout, PARAM, CONFIG, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
    [psi, CCR, CCZ, DELGE, RCCS, ZCCS, FL2, BZ2, BR2] = INTER(PARAM, 0, FFout, ExtCOIL,...
    SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % 結果を表示
    if CONFIG.ShowFig
        show_results(PARAM, CONFIG, CCR, CCZ, CCSDAT, REF, psi,...
            SENSOR_TPRB, SENSOR_FLXLP, WALL,...
            ExtCOIL, DISF, t, Ip)
    end
end
