clear all

aaa = load('vars_result_良い', "FFout");
FFoutGood = aaa.FFout;

type = "sol";
[PARAM, CONFIG] = define_params(type);
[FL, BZ, BR, ExtCOIL, REF, JEDDY, POS] = load_num_sol(PARAM, CONFIG); % OK
CONFIG.ShowFig = 0;
CONFIG.ShowFig2 = 1;

% 3つの入力データで調べる→手動

gosa = [0, 0.03];
jigen = [0, 1];

fji = figure();
cnt = 0;

% 誤差なし一致なし
% 誤差なし一致あり
% 誤差あり一致なし
% 誤差あり一致あり

% 誤差の有無
for G = 1:2

    PARAM.SIGM = gosa(G);

    % 次元の一致の有無
    for J = 1:2

        cnt = cnt + 1;
        CONFIG.DevideFlux = jigen(J);

        % 各パラメータを設定
        [FL, BZ, BR, ExtCOIL, REF, JEDDY, POS] = load_num_sol(PARAM, CONFIG); % OK
        [PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP] = make_sensor(PARAM, CONFIG, BZ, BR, FL, POS);
        CCSDAT = make_CCS(PARAM);
        FF = make_FF(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT);
        % WALL = load_wall();
        WALL = loadwalldata2(PARAM, CONFIG);
        % 順問題の解行列FF、CCS点と各センサー及び自分以外のCCS点との関係式行列AAを作成
        [FC, BR, BZ, PSIFLX, PSIC, AA, FF] = FORM(PARAM, CONFIG, FF, ...
            ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);
        % 重み付け
        [FF, FC, AA, factors] = weight(CONFIG, FF,FC,AA,SENSOR_NPRB,SENSOR_TPRB,SENSOR_FLXLP,CCSDAT);
        % 逆問題を解く
        FFout = tsvd(PARAM, CONFIG, AA, FF, FC, ...
        SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);
        % 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
        [psi, CCR, CCZ, DELGE, RCCS, ZCCS, FL2, BZ2, BR2] = INTER(PARAM, 0, FFout, ExtCOIL,...
        SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

        err_LCFS(cnt) = evaluate_LCFS(psi, REF, PARAM, CONFIG, CCR, CCZ, 0);
        residual(cnt) = norm(AA * FFout - FF');

        figure(fji);
        subplot(1,2,G)
        if J == 1
            hold on
            v = linspace(-20, 20, 101);
            REFmini = imresize(REF.Flux, size(psi));
            contour(CCR, CCZ, REFmini*1000, v)
            contour(CCR, CCZ, REFmini*1000, [0 0], 'LineColor', 'k', 'LineWidth', 2);
            contour(CCR, CCZ, psi, [0 0], 'LineColor', 'c', 'LineWidth', 2);
            plot(WALL.RWALL_list, WALL.ZWALL_list, '-k'); % 容器壁 VacuumVesselMeshPoints
            xlabel({'r (m)'});
            ylabel({'z (m)'});
            % title("Reconstructed flux")
            axis equal
            hold off
        else
            hold on
            contour(CCR, CCZ, psi, [0 0], 'LineColor', 'm', 'LineWidth', 2);
            hold off
        end
        figure()
        contour(CCR, CCZ, psi);
    end
end

disp(err_LCFS)
disp(residual)