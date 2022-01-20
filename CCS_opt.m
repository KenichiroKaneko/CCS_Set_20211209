close all
clear all
% 入力となるクソデカセンサー配置群を作成
type = 'sol';

% [PARAM, CONFIG] = define_params(type);
% PARAM.MP_pos = "/CCS_sensor_position_opt.txt";
% PARAM.FL_pos = "/CCS_sensor_position_opt.txt";
% PARAM.KUP = 0;
% CONFIG.RevSenPos = 1;
% CONFIG.ShowFig = 1;
% [FL, BZ, BR, ExtCOIL, REF, JEDDY, POS] = load_num_sol(PARAM, CONFIG); % OK
% t = 0; Ip = 0;
% [PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP] = make_sensor(PARAM, CONFIG, BZ, BR, FL, POS);
% CCSDAT = make_CCS(PARAM);
% FF = make_FF(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT);
% WALL = loadwalldata2(PARAM, CONFIG);
% % 順問題の解行列FF、CCS点と各センサー及び自分以外のCCS点との関係式行列AAを作成
% [FC, BR, BZ, PSIFLX, PSIC, AA, FF] = FORM(PARAM, CONFIG, FF, ...
% ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);
% % 重み付け
% [FF, FC, AA, factors] = weight(FF,FC,AA,SENSOR_NPRB,SENSOR_TPRB,SENSOR_FLXLP,CCSDAT);
% % 逆問題を解く
% FFout = tsvd(PARAM, CONFIG, AA, FF, FC, ...
%     SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);
% save('vars_verybig2033')

CONFIG.RevSenPos = 0;
CONFIG.ShowFig = 0;
PARAM.KUP = 0;

GA = opt_GA();
GA.evolve();

figure()
plot(GA.err)

for i = 1:GA.N
    errs(i) = GA.pool(i).err;
end

[M I] = min(errs);
Gene = GA.pool(I).gene;

BZ_i = find(Gene(1 : GA.pool(I).BZ_gene_num) == 1);
BZ_t = sprintf(repmat('% 3d', [1, length(BZ_i)]), BZ_i);
FL_i = find(Gene(GA.pool(I).BZ_gene_num+1 : end) == 1);
FL_t = sprintf(repmat('% 3d', [1, length(FL_i)]), FL_i);
err_t = sprintf('   %.5e', GA.pool(I).err);
disp([err_t BZ_t FL_t]);

state = '1月18日の深夜にスタートした50個の遺伝子を4000代までやったやつ';

save('vars_GA_0118')
% error('error description', A1)


