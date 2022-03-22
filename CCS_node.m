% 渦電流ノード点の最適な位置を探索

close all; clear all;

tStart = tic;
f1 = figure();
j = 0;

for i = 0:3^7-1
    % 3進数を使って渦電流ノード点を全探索
    str = dec2base(i,3,7);
    node_list = str2num(str');
    node_list = node_list' + 1;

    node_num = sum(3 * [node_list node_list(1:end-1)]);
    if node_num > 57
        continue
    end
    j = j+1;
    
    % 渦電流
    [PARAM, CONFIG] = define_params('sol');
    PARAM.eddyNodes = node_list;
    [err_LCFS, resid] = CCS_PARAM(PARAM, CONFIG);
    err_LCFS_list(j) = err_LCFS;
    residual_list(j) = resid;

    figure(f1)
    hold on
    plot(j, err_LCFS, 'mo')
    plot(j, resid, 'co')
    hold off

    disp([j node_list])
    if rem(j,10) == 0
        disp(toc(tStart))
    end
end

figure()
hold on
plot(1:j, err_LCFS_list)
[m, I] = min(err_LCFS_list);
plot(I, err_LCFS_list(I), 'o')
disp(m)
disp(I)
disp(dec2base(I, 2, 7))
plot(1:j, residual_list)
[m, I] = min(residual_list);
plot(I, residual_list(I), 'o')
disp(['経過時間', num2str(toc(tStart))])

