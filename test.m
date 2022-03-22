% 新たなLCFSの評価
clear all;
close all;
load('vars_result')

evaluate_LCFS(psi, REF, PARAM, CONFIG, CCR, CCZ, 0)
% 再構成結果の表示
% psiBig = imresize(psi, size(REF.Flux));
v = linspace(-20, 20, 101);


REF.Flux = imresize(REF.Flux, size(psi));
REF.Flux(:,1:2) = 1;
REF01 = REF.Flux<0;
se = strel('disk',2);
% クロージング処理
REF01 = imclose(REF01, se);
% 収縮処理
REF01 = imerode(REF01, se);
REF.Flux(REF01) = -1;

figure()
imshow(REF01)

REF.Flux(REF01) = 0;
figure()
contour(CCR, CCZ, REF.Flux*1000, v, 'm')


error('error description', A1)

subplot(1, 2, 2);
hold on

contour(REF.R, REF.Z, REF.Flux*1000, v, 'm'); % 正解;
contour (REF.R, REF.Z, REF.Flux, [0 0], 'm', 'LineWidth', 2);
contour(REF.R, REF.Z, psiBig*1000, v, 'k');
contour(REF.R, REF.Z, psiBig, [0 0], 'c', 'LineWidth', 2)

plot(CCSDAT.RGI, CCSDAT.ZGI, CCSDAT.RGI, -CCSDAT.ZGI, 'color','y','LineWidth', 2);
plot(WALL.RWALL_list, WALL.ZWALL_list, '-k'); % 容器壁 VacuumVesselMeshPoints
hold off
xlabel({'r (m)'});
ylabel({'z (m)'});
title("Reconstructed flux")
axis equal

CONFIG.ShowFig = 1;

% err_LCFS = evaluate_LCFS(psi, REF, PARAM, CONFIG, CCR, CCZ,0)

psi_add = (REF.Flux - psiBig);

figure()
hold on
psiBig = psiBig + psi_add;
contour(REF.R, REF.Z, psiBig*1000, v, 'k');
contour(REF.R, REF.Z, psiBig, [0 0], 'c', 'LineWidth', 2)
contour(REF.R, REF.Z, REF.Flux*1000, v, 'm--'); % 正解;
contour (REF.R, REF.Z, REF.Flux, [0 0], 'm', 'LineWidth', 2);
plot(CCSDAT.RGI, CCSDAT.ZGI, CCSDAT.RGI, -CCSDAT.ZGI, 'color','y','LineWidth', 2);
plot(WALL.RWALL_list, WALL.ZWALL_list, '-k'); % 容器壁 VacuumVesselMeshPoints
hold off
xlabel({'r (m)'});
ylabel({'z (m)'});
title("Reconstructed flux")
axis equal

save('psi_add', 'psi_add', 'CCR', 'CCZ')