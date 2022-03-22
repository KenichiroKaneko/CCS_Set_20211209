% result_2033_4 residual, err_FFoutが正しい
% result_2033_1,2 FFとか間違えたからダメデータ
% result_2033_3％誤差 3％誤差

% result_2033_err1 誤差の配列をしまってある

% めも
% 一回目をほどほどな回(誤差なし)
% 2回目を打切り項数を動かす

close all; clear all;
type = 'exp';
type = 'sol';

REF = load('vars_result_良い')
REF1.Flux = REF.psi;
REF1.R = REF.CCR;
REF1.Z = REF.CCZ;
FF1 = REF.FF;
FFout1 = REF.FFout;
AA1 = REF.AA;
FC1 = REF.FC;
PARAM = REF.PARAM;
CONFIG = REF.CONFIG;
SENSOR_TPRB = REF.SENSOR_TPRB;
SENSOR_NPRB = REF.SENSOR_NPRB;
SENSOR_FLXLP = REF.SENSOR_FLXLP;
ExtCOIL = REF.ExtCOIL;
CCSDAT = REF.CCSDAT;
WALL = REF.WALL;

CONFIG.ShowFig = 0;
CONFIG.ShowFig2 = 0;

tic
for I = 1:25
    I
    PARAM.KUP = I + 50;
    AA2 = AA1;
    FF2 = (AA1 * FFout1)';
    FC2 = FC1;

    % 逆問題を解く
    FFout2 = tsvd(PARAM, CONFIG, AA2, FF2, FC2, ...
        SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
    [psi2, CCR, CCZ, DELGE, RCCS, ZCCS, FL2, BZ2, BR2] = INTER(PARAM, 0, FFout2, ExtCOIL,...
    SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % LCFS残差1段目の良好な再構成結果と2段目との比較
    err_LCFS(I) = evaluate_LCFS(psi2, REF1, PARAM, CONFIG, CCR, CCZ, 0);
    % FF-AA*X 残差
    residual(I) = norm(AA2 * FFout2 - FF2);
    % ccs,node残差???
    err_FFout(I) = norm(FFout2 - FFout1);
    
    save("./cal/result_2033_良好な1段目との比較/"+num2str(I), "psi2", "FF2", "FFout2", "AA2");
    % save("./cal/result_2033_err1/"+num2str(I), "FF", "FFout", "AA", "FFout2");

end
toc

save("./cal/result_2033_良好な1段目との比較/errs")

