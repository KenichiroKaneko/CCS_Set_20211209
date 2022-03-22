function CCS_UTST()
    %
    close all; 
    clear all;
    type = 'exp';
    % type = 'sol';
    [PARAM, CONFIG] = define_params(type);
    CONFIG.ShowFig = 1;
    CONFIG.ShowFig2 = 1;

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

    if CONFIG.WithExpCoil == 1
        [FL_, BZ_, ExtCOIL, Ip_] = load_lizzie(PARAM, CONFIG);
        t = round((str2double(PARAM.time_CCS)) / 0.5);
        ExtCOIL.I(3:10) = ExtCOIL.I_sig(1:8, t);
    end

    save([PARAM.output_file_directory PARAM.num_sol_name],...
    'CONFIG', 'FL', 'BZ', 'BR', 'ExtCOIL', 'REF', 'JEDDY', 'POS');
    % save('ExtCOIL', 'ExtCOIL');
    % figure()
    % subplot(2,1,1)
    % hold on
    % plot(SENSOR_FLXLP.FLXLP, 'o-')
    % plot(SENSOR_FLXLP.FLXLP ./ SENSOR_FLXLP.R, 'o-')
    % legend('次元の一致なし','次元の一致あり')
    % subplot(2,1,2)
    % plot(SENSOR_TPRB.TPRB,'o-')
    % disp("分散")
    % disp("全信号   " + num2str(var([SENSOR_TPRB.TPRB SENSOR_FLXLP.FLXLP])))
    % disp("磁場内側 " + num2str(var(SENSOR_TPRB.TPRB(1:18))))
    % disp("磁場外側 " + num2str(var(SENSOR_TPRB.TPRB(20:39))))
    % disp("磁束内側 " + num2str(var(SENSOR_FLXLP.FLXLP(1:19))))
    % disp("磁束外側 " + num2str(var(SENSOR_FLXLP.FLXLP(20:34))))
    % disp("磁場     " + num2str(var(SENSOR_TPRB.TPRB)))
    % disp("磁束     " + num2str(var(SENSOR_FLXLP.FLXLP)))
    % disp("平均")
    % disp("全信号   " + num2str(mean([SENSOR_TPRB.TPRB SENSOR_FLXLP.FLXLP])))
    % disp("磁場内側 " + num2str(mean(SENSOR_TPRB.TPRB(1:18))))
    % disp("磁場外側 " + num2str(mean(SENSOR_TPRB.TPRB(20:39))))
    % disp("磁束内側 " + num2str(mean(SENSOR_FLXLP.FLXLP(1:19))))
    % disp("磁束外側 " + num2str(mean(SENSOR_FLXLP.FLXLP(20:34))))
    % disp("磁場     " + num2str(mean(SENSOR_TPRB.TPRB)))

    CCSDAT = make_CCS(PARAM);
    FF = make_FF(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT);
    % WALL = load_wall();
    WALL = loadwalldata2(PARAM, CONFIG);

    figure()
    hold on
    plot(WALL.REV, WALL.ZEV, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8) % VacuumVesselNodePoints
    plot(WALL.RSEC, WALL.ZSEC, 'mo', 'MarkerSize', 16) % VacuumVesselSegments
    plot(WALL.RWALL_list, WALL.ZWALL_list, '-k'); % 容器壁 VacuumVesselMeshPoints
    % legend();
    xlabel({'r (m)'});
    ylabel({'z (m)'});
    % % xlim([0, 0.7])
    % % ylim([-1, 1])
    % axis equal

    
    % 順問題の解行列FF、CCS点と各センサー及び自分以外のCCS点との関係式行列AAを作成
    [FC, BR, BZ, PSIFLX, PSIC, AA, FF] = FORM(PARAM, CONFIG, FF, ...
        ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % 重み付け
    [FF, FC, AA, factors] = weight(CONFIG, FF,FC,AA,SENSOR_NPRB,SENSOR_TPRB,SENSOR_FLXLP,CCSDAT);
    % 逆問題を解く
    FFout = tsvd(PARAM, CONFIG, AA, FF, FC, ...
        SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % % プロット用に渦電流を計算する
    DISF = EDDYP(FFout, PARAM, CONFIG, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
    [psi, CCR, CCZ, DELGE, RCCS, ZCCS, FL2, BZ2, BR2] = INTER(PARAM, 0, FFout, ExtCOIL,...
        SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    save('vars_result_良い')
    save('vars_result')

    % 結果を表示
    if CONFIG.ShowFig
        show_results(PARAM, CONFIG, CCR, CCZ, CCSDAT, REF, psi,...
            SENSOR_TPRB, SENSOR_FLXLP, WALL,...
            ExtCOIL, DISF, t, Ip, FFout)
    end

    err_LCFS = evaluate_LCFS(psi, REF, PARAM, CONFIG, CCR, CCZ, 0);
    txt = sprintf('%.3e', err_LCFS);
    disp(['LCFS残差: ' txt])
    
end
