close all; clear all;
type = 'sol';
[PARAM, CONFIG] = define_params(type);
[FL, BZ, BR, ExtCOIL, REF, JEDDY, POS] = load_num_sol(PARAM, CONFIG); % OK
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
[FF, FC, AA, factors] = weight(CONFIG, FF,FC,AA,SENSOR_NPRB,SENSOR_TPRB,SENSOR_FLXLP,CCSDAT);

save('vars_inputs');

figure('Name', 'LCFS残差')
ani = animatedline;

for i = 1:26
    break
    KUP = i + 49;
    % load('vars_inputs');
    % PARAM.KUP = KUP;
    % % 逆問題を解く
    % FFout = tsvd(PARAM, CONFIG, AA, FF, FC, ...
    % SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);
    % % 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
    % [psi, CCR, CCZ, DELGE, RCCS, ZCCS, FL2, BZ2, BR2] = INTER(PARAM, 0, FFout, ExtCOIL,...
    %     SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % load("cal/result_2033_良好なKUPを探索誤差なし" + num2str(i))
    load("cal/result_2033_良好な1段目との比較/"+num2str(i))
    psiBig = imresize(psi, size(REF.Flux));
    load('psi_add', "psi_add", "CCR", "CCZ")
    psi = psiBig + psi_add;
    
    % LCFS残差
    err_LCFS(i) = evaluate_LCFS(psi, REF, PARAM, CONFIG, CCR, CCZ, 0);
    % FF-AA*X 残差
    residual(i) = norm(AA * FFout - FF);

    % save("cal/result_2033_良好なKUPを探索誤差なし" + num2str(i))

    addpoints(ani, i, err_LCFS(i))
    drawnow
    disp(KUP);
end
close all; clear all;
w = ["max", "norm", "10", "100", "no"];

type = 'sol';
[PARAM, CONFIG] = define_params(type);
CONFIG.ShowFig = 0;
CONFIG.ShowFig2 = 0;
f1 = figure();
hold on
hold off

for i = 1:5
    CONFIG.Nmrz = w(i); 

    [FL, BZ, BR, ExtCOIL, REF, JEDDY, POS] = load_num_sol(PARAM, CONFIG); % OK
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
    [FF, FC, AA, factors] = weight(CONFIG, FF, FC, AA, SENSOR_NPRB, SENSOR_TPRB, SENSOR_FLXLP, CCSDAT);

    % 逆問題を解く
    FFout = tsvd(PARAM, CONFIG, AA, FF, FC, ...
    SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % プロット用に渦電流を計算する
    DISF = EDDYP(FFout, PARAM, CONFIG, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);
    % 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
    [psi, CCR, CCZ, DELGE, RCCS, ZCCS, FL2, BZ2, BR2] = INTER(PARAM, 0, FFout, ExtCOIL,...
    SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % LCFS残差
    err_LCFS(i) = evaluate_LCFS(psi, REF, PARAM, CONFIG, CCR, CCZ, 0);
    figure(f1)
    hold on
    subplot(1,5,i)
    v = linspace(-20, 20, 101);
    contour(REF.R, REF.Z, REF.Flux*1000, v, 'k'); % 正解;
    contour(REF.R, REF.Z, REF.Flux, [0 0], 'LineColor', 'k', 'LineWidth', 2);
    contour(CCR, CCZ, psi, [0 0], 'LineColor', 'm', 'LineWidth', 2);
    plot(CCSDAT.RGI, CCSDAT.ZGI, CCSDAT.RGI, -CCSDAT.ZGI, 'color','y','LineWidth', 2);
    plot(WALL.RWALL_list, WALL.ZWALL_list, '-k'); % 容器壁 VacuumVesselMeshPoints
    xlabel({'r (m)'});
    ylabel({'z (m)'});
    % title("Reconstructed flux")
    axis equal
    hold off

end




% figure()
% hold on
% xlabel("LCFS誤差")
% ylabel("残差 ||Dp^* - g||")
% plot(err_LCFS, residual, 'o')
% plot(err_LCFS(20), residual(20), '*')


% for i = 1:length(err_LCFS)
%     text(err_LCFS(i), (residual(i)), ['\leftarrow', num2str(i + 49)])
% end

% figure()
% hold on
% [m, I] = min(residual)
% plot(residual)
% plot(I, residual(I), 'o')

% figure()
% plot(err_LCFS)
