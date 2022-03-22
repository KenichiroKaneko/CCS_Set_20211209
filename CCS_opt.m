close all
clear all
% 入力となるクソデカセンサー配置群を作成
type = 'sol';

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
Gene = GA.goodgene;

BZ_i = find(Gene(1 : GA.pool(I).BZ_gene_num) == 1);
BZ_t = sprintf(repmat('% 3d', [1, length(BZ_i)]), BZ_i);
FL_i = find(Gene(GA.pool(I).BZ_gene_num+1 : end) == 1);
FL_t = sprintf(repmat('% 3d', [1, length(FL_i)]), FL_i);
err_t = sprintf('   %.5e', GA.pool(I).err);
disp([err_t BZ_t FL_t]);

state = '3月1日三井GH';

save('vars_GA_0301')
% error('error description', A1)


