% CCSセンサーでどのセンサーがどれだけ重要か知る
close all;
clear all;

DISF = 0;
load('vars_exp')

% 元になるやつ
AA_ = AA;
FF_ = FF;
FFout_ = FFout;
SENSOR_NPRB_ = SENSOR_NPRB;
SENSOR_TPRB_ = SENSOR_TPRB;
SENSOR_FLXLP_ = SENSOR_FLXLP;
NCCN = CCSDAT.NCCN;
REF.Z = CCZ;
REF.R = CCR;
REF.Flux = psi;

CONFIG.ShowFig = 0;

for i = 1:SENSOR_TPRB_.NUM + SENSOR_FLXLP_.NUM
% for i = 1:1
    % 解からもう一度逆問題を解く
    FF_ = AA_ * FFout_;
    FF_ = FF_';
    if i <= SENSOR_TPRB_.NUM
        PARAM.dead_BZ = [i];
        PARAM.dead_FL = [];
    else
        PARAM.dead_BZ = [];
        PARAM.dead_FL = [i - SENSOR_TPRB_.NUM]; 
    end
    [FF, AA, SENSOR_TPRB, SENSOR_FLXLP] = remake_sensor(PARAM, FF_, AA_, SENSOR_FLXLP_, SENSOR_TPRB_);

    % 重み付け
    [FF, FC, AA, factors] = weight(FF, FC, AA, SENSOR_NPRB, SENSOR_TPRB, SENSOR_FLXLP, CCSDAT);

    % 逆問題を解く
    FFout = tsvd(PARAM, CONFIG, AA, FF, FC, ...
        SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
    [psi, CCR, CCZ, DELGE, RCCS, ZCCS, FL2, BZ2, BR2] = INTER(PARAM, 0, FFout, ExtCOIL,...
    SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % show_results(PARAM, CONFIG, CCR, CCZ, CCSDAT, REF, psi, SENSOR_TPRB, SENSOR_FLXLP, WALL,...
    % ExtCOIL, DISF, t, Ip)

    err_LCFS(i) = evaluate_LCFS(psi, REF, PARAM, CCR, CCZ, i);
    err_residual(i, :) = abs(FF - (AA * FFout)');
    err_node(i, :) = abs(FFout - FFout_);
    FFout_list(i, :) = FFout;
end

figure('Name', 'LCFSの誤差')
subplot(1,2,1);
plot(err_LCFS(1:SENSOR_TPRB.NUM));title('BZ');
subplot(1,2,2);
plot(err_LCFS(SENSOR_TPRB.NUM+1:end));title('FL');

figure('Name', 'residual')
subplot(1,2,1); hold on;
plot(sum(err_residual(1:SENSOR_TPRB.NUM, :), 2));title('BZ residual')
subplot(1,2,2); hold on;
plot(sum(err_residual(SENSOR_TPRB.NUM+1:end, :), 2));title('FL residual')

figure('Name', 'node')
subplot(2,2,1);
plot(sum(err_node(1:SENSOR_TPRB.NUM, 1:NCCN), 2));title('BZ CCS');
subplot(2,2,2);
plot(sum(err_node(1:SENSOR_TPRB.NUM, NCCN+1:end), 2));title('BZ eddy');
subplot(2,2,3);
plot(sum(err_node(SENSOR_TPRB.NUM+1:end, 1:NCCN), 2));title('FL CCS');
subplot(2,2,4);
plot(sum(err_node(SENSOR_TPRB.NUM+1:end, NCCN+1:end), 2));title('BZ eddy');

save('vars_roop_合体後1')