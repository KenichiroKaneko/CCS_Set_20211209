function [err_LCFS residual] = CCS_PARAM(PARAM, CONFIG)
    % [FL, BZ, BR, ExtCOIL, REF, JEDDY, POS] = load_num_sol(PARAM, CONFIG); % OK
    load([PARAM.output_file_directory PARAM.num_sol_name]);
    CONFIG.ShowFig = 0;
    CONFIG.ShowFig2 = 0;
    % 各パラメータを設定
    [PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP] = make_sensor(PARAM, CONFIG, BZ, BR, FL, POS);
    CCSDAT = make_CCS(PARAM);
    WALL = loadwalldata2(PARAM, CONFIG);
    FF = make_FF(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT);
    % 順問題の解行列FF、CCS点と各センサー及び自分以外のCCS点との関係式行列AAを作成
    [FC, BR, BZ, PSIFLX, PSIC, AA, FF] = FORM(PARAM, CONFIG, FF, ...
        ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);
    % 重み付け
    [FF, FC, AA, factors] = weight(CONFIG, FF,FC,AA,SENSOR_NPRB,SENSOR_TPRB,SENSOR_FLXLP,CCSDAT);
    % 逆問題を解く
    FFout = tsvd(PARAM, CONFIG, AA, FF, FC, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);
    % 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
    [psi, CCR, CCZ, DELGE, RCCS, ZCCS, FL2, BZ2, BR2] = INTER(PARAM, 0, FFout, ExtCOIL,...
        SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % REFmini = imresize(REF.Flux, size(psi));
    % v = linspace(-20, 20, 101);
    % figure()
    % subplot(1,2,1)
    % hold on
    % contour(CCR, CCZ, REFmini*1000, v)
    % contour(CCR, CCZ, REFmini*1000, [0 0], 'LineColor', 'm', 'LineWidth', 2);
    % contour(CCR, CCZ, psi, [0 0], 'k', 'LineWidth', 2)
    % plot(WALL.REV, WALL.ZEV, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8) % VacuumVesselNodePoints
    % plot(WALL.RSEC, WALL.ZSEC, 'mo', 'MarkerSize', 16) % VacuumVesselSegments
    % plot(WALL.RWALL_list, WALL.ZWALL_list, '-k'); % 容器壁 VacuumVesselMeshPoints
    % % legend();
    % xlabel({'r (m)'});
    % ylabel({'z (m)'});
    % axis equal
    % hold off;
    % subplot(1,2,2)
    % hold on
    % contour(CCR, CCZ, REFmini*1000, v, 'm'); % 正解;
    % contour(CCR, CCZ, REFmini*1000, [0 0], 'm', 'LineWidth', 2);
    % contour(CCR, CCZ, psi*1000, v,'k');
    % contour(CCR, CCZ, psi, [0 0], 'LineColor', 'c', 'LineWidth', 2);
    % plot(WALL.RWALL_list, WALL.ZWALL_list, '-k'); % 容器壁 VacuumVesselMeshPoints
    % hold off
    % xlabel({'r (m)'});
    % ylabel({'z (m)'});
    % title("Reconstructed flux")
    % axis equal
    residual = norm(AA * FFout - FF');
    err_LCFS = evaluate_LCFS(psi, REF, PARAM, CONFIG, CCR, CCZ, 0);
end