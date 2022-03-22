% 3つのデータで調べる
% 打切り項数を動かす
f1 = figure()
for K = 1:75
    disp(K)
    j = K;
    [PARAM, CONFIG] = define_params('sol');
    PARAM.KUP = K;
    [err_LCFS, resid] = CCS_PARAM(PARAM, CONFIG);
    err_LCFS_list(j) = err_LCFS;
    residual_list(j) = resid;

    figure(f1)
    hold on
    plot(j, err_LCFS, 'mo')
    plot(j, resid, 'co')
    hold off
    % [PARAM, CONFIG] = define_params(type);
    % [FL, BZ, BR, ExtCOIL, REF, JEDDY, POS] = load_num_sol(PARAM, CONFIG); % OK
    % CONFIG.ShowFig = 0;
    % CONFIG.ShowFig2 = 0;
    % PARAM.KUP = K + 50;
    % % 各パラメータを設定
    % [PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP] = make_sensor(PARAM, CONFIG, BZ, BR, FL, POS);
    % CCSDAT = make_CCS(PARAM);
    % FF = make_FF(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT);
    % WALL = loadwalldata2(PARAM, CONFIG);
    % % 順問題の解行列FF、CCS点と各センサー及び自分以外のCCS点との関係式行列AAを作成
    % [FC, BR, BZ, PSIFLX, PSIC, AA, FF] = FORM(PARAM, CONFIG, FF, ...
    %     ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);
    % % 重み付け
    % [FF, FC, AA, factors] = weight(CONFIG, FF,FC,AA,SENSOR_NPRB,SENSOR_TPRB,SENSOR_FLXLP,CCSDAT);
    % % 逆問題を解く
    % FFout = tsvd(PARAM, CONFIG, AA, FF, FC, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);
    % % 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
    % [psi, CCR, CCZ, DELGE, RCCS, ZCCS, FL2, BZ2, BR2] = INTER(PARAM, 0, FFout, ExtCOIL,...
    %     SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % save("cal/result_1000_1段目でLを検証2022年2月28日/"+num2str(K), 'FF', 'AA', 'FC', 'FFout', 'psi', 'CCR', 'CCZ', 'PARAM', 'WALL', 'CCSDAT', 'CONFIG')
    
    % err_LCFS(K) = evaluate_LCFS(psi, REF, PARAM, CONFIG, CCR, CCZ, 0);
end

save("cal/result_2033_LcurveKUP")
figure()
plot(51:75, err_LCFS)
hold on 
[m I] = min(err_LCFS)
plot(I+50, m, 'o')
