
close all; clear all;
type = 'exp';
type = 'sol';
[PARAM, CONFIG] = define_params(type);
% 各パラメータを設定
[FL, BZ, BR, ExtCOIL, REF, JEDDY, POS] = load_num_sol(PARAM, CONFIG); 
[PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP] = make_sensor(PARAM, CONFIG, BZ, BR, FL, POS);
CCSDAT = make_CCS(PARAM);
FF = make_FF(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT);
WALL = loadwalldata2(PARAM, CONFIG);
% 順問題の解行列FF、CCS点と各センサー及び自分以外のCCS点との関係式行列AAを作成
[FC, BR, BZ, PSIFLX, PSIC, AA, FF] = FORM(PARAM, CONFIG, FF, ...
    ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);
% 重み付け
[FF, FC, AA, factors] = weight(CONFIG, FF,FC,AA,SENSOR_NPRB,SENSOR_TPRB,SENSOR_FLXLP,CCSDAT);

% FFout75 = load('./cal/KUP75', 'FFout');

CONFIG.ShowFig = 0;
CONFIG.ShowFig2 = 0;

tic
for I = 1:75
    PARAM.KUP = I;

    % 逆問題を解く
    FFout = tsvd(PARAM, CONFIG, AA, FF, FC, ...
        SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);
    % 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
    [psi, CCR, CCZ, DELGE, RCCS, ZCCS, FL2, BZ2, BR2] = INTER(PARAM, 0, FFout, ExtCOIL,...
    SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % LCFS残差
    err_LCFS(I) = evaluate_LCFS(psi, REF, PARAM, CONFIG, CCR, CCZ, 0);
    % FF-AA*X 残差
    residual(I) = norm(AA * FFout - FF');

    % 解からもう一度逆問題を解く
    FF = (AA * FFout)';
    PARAM.dead_FL = []; PARAM.dead_BZ = [];
    % 重み付け
    [FF, FC, AA, factors] = weight(CONFIG,FF,FC,AA,SENSOR_NPRB,SENSOR_TPRB,SENSOR_FLXLP,CCSDAT);
    % 逆問題を解く
    FFout2 = tsvd(PARAM, CONFIG, AA, FF, FC, ...
        SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    ccs,node残差???
    err_FFout(I) = norm(FFout2 - FFout);
    
    save("./cal/result_2033/"+num2str(I), "psi", "FF", "FFout", "AA");
    
end
toc

err_LCFS
residual

save("./cal/result_2033/errs")

load("./cal/result_2033/errs")
figure()
plot(err_LCFS)


