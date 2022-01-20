% CCS
clear all; close all;
type = 'sol';
[PARAM, CONFIG] = define_params(type);
PARAM.MP_pos = "/CCS_sensor_position_opt.txt";
PARAM.FL_pos = "/CCS_sensor_position_opt.txt";
PARAM.KUP = 0;
CONFIG.RevSenPos = 1;
CONFIG.ShowFig = 1;
[FL, BZ, BR, ExtCOIL, REF, JEDDY, POS] = load_num_sol(PARAM, CONFIG); % OK
t = 0; Ip = 0;
[PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP] = make_sensor(PARAM, CONFIG, BZ, BR, FL, POS);
CCSDAT = make_CCS(PARAM);
WALL = loadwalldata2(PARAM, CONFIG);

FLgene = zeros(1, 52);
BZgene = zeros(1, 52);

% ind_BZ = [2  4  6  7 13 16 17 18 28 30 31 37 39 40 43 47 49 52  ];
% ind_FL = [2  4 12 16 20 21 22 24 31 34 36 37 38 40 45 46 49];
% ind_BZ = [1  2  4  5  7  8 13 14 22 28 33 36 45 46 47 48 50 51];
% ind_FL = [3  5  6  8 11 16 18 20 26 27 29 35 39 43 44 46 49 50 51];
ind_BZ = [1  2  4  5  8  9 11 13 14 16 19 22 23 28 41 42 45 48 49];
ind_FL = [3  4  6 11 12 16 19 26 29 32 35 36 37 42 44 45 46 47 49];
% ind_BZ = [2 3:3:51 52];
% ind_FL = [2 3:3:51 52];
FLgene(ind_FL) = 1;
BZgene(ind_BZ) = 1;
ind_BZ = [1, BZgene, 1, fliplr(BZgene)];
ind_FL = [1, FLgene, 1, fliplr(FLgene)];
ind_BZ = find(ind_BZ == 1);
ind_FL = find(ind_FL == 1);
% SENSORを作り直す
% PARAM.KUP = 65;
% SENSOR_TPRB.R = SENSOR_TPRB.R(ind_BZ);
% SENSOR_TPRB.Z = SENSOR_TPRB.Z(ind_BZ);
% SENSOR_TPRB.TET = SENSOR_TPRB.TET(ind_BZ);
% SENSOR_TPRB.TPRB = SENSOR_TPRB.TPRB(ind_BZ);
% SENSOR_TPRB.NUM = length(SENSOR_TPRB.TPRB);
% SENSOR_FLXLP.R = SENSOR_FLXLP.R(ind_FL);
% SENSOR_FLXLP.Z = SENSOR_FLXLP.Z(ind_FL);
% SENSOR_FLXLP.TET = SENSOR_FLXLP.TET(ind_FL);
% SENSOR_FLXLP.FLXLP = SENSOR_FLXLP.FLXLP(ind_FL);
% SENSOR_FLXLP.NUM = length(SENSOR_FLXLP.FLXLP);


% 順問題の解行列FF、CCS点と各センサー及び自分以外のCCS点との関係式行列AAを作成
FF = make_FF(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT);
[FC, BR, BZ, PSIFLX, PSIC, AA, FF] = FORM(PARAM, CONFIG, FF, ...
ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

% 重み付け
[FF, FC, AA, factors] = weight(FF,FC,AA,SENSOR_NPRB,SENSOR_TPRB,SENSOR_FLXLP,CCSDAT);

% 逆問題を解く
FFout = tsvd(PARAM, CONFIG, AA, FF, FC, ...
SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);
err = sum(abs(FF' - (AA * FFout)));

% プロット用に渦電流を計算する
DISF = EDDYP(FFout, PARAM, CONFIG, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

% % 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
[psi, CCR, CCZ, DELGE, RCCS, ZCCS, FL2, BZ2, BR2] = INTER(PARAM, 0, FFout, ExtCOIL,...
SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

show_results(PARAM, CONFIG, CCR, CCZ, CCSDAT, REF, psi,...
            SENSOR_TPRB, SENSOR_FLXLP, WALL,...
            ExtCOIL, DISF, t, Ip, FFout)


% % 誤差評価する
% err = evaluate_LCFS(psi, REF, PARAM, CONFIG, CCR, CCZ, 0)

% load('vars_GA_0118')
% figure()
% plot(GA.err)
% xlabel('世代')
% ylabel('適応度（小さいほど良い）')
