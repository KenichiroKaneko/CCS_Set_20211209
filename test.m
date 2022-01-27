% 新たなLCFSの評価
clear all;
close all;
load('vars_result_KUP75')

CONFIG.ShowFig = 1;
err_LCFS = evaluate_LCFS(psi, REF, PARAM, CONFIG, CCR, CCZ,0)

% % return
% CONFIG.ShowFig = 1;
% MSE = EVALUATE(psi, REF, PARAM, CONFIG, CCSDAT, CCR, CCZ)
% % err_LCFS = evaluate_LCFS(psi, REF, PARAM, CONFIG, CCR, CCZ, 0)
% error('error description', A1)
% l = 0.6;
% s = 1000;
% ignore = round(s * 0.1);
% n = 32;
% zhalf = round(length(REF.Z) / 2);
% figure()
% hold on
% for c = 1:2
%     if c == 1
%         range = 1:zhalf;
%         zplus = 0;
%     else
%         range = zhalf+1:length(REF.Z);
%         zplus = zhalf;
%     end
%     % 磁気軸
%     [m r0] = min(min(REF.Flux(range, :)));
%     [m z0] = min(REF.Flux(range, r0));
%     R0(c) = REF.R(r0);
%     Z0(c) = REF.Z(z0 + zplus);
%     plot(R0(c), Z0(c), 'o')
%     for i = 1:n
%         theta = 2 * pi / n * (i - 1) + pi / n;
%         r = l * sin(theta); z = l * cos(theta);
%         line_r = [R0(c), R0(c) + r]; line_z = [Z0(c), Z0(c) + z];
%         % ある直線が作るpsiの断面を取得 
%         [cx_ref, cy_ref, c_ref] = improfile(REF.R, REF.Z(range), REF.Flux(range, :), line_r, line_z, s, 'bilinear');
%         c_ref(1:ignore) = -0.001;
%         % LCFS(psi=0)との交点を探索
%         [m, J] = min(abs(c_ref));
%         I_ref(c, i) = J;
%         plot3(cx_ref, cy_ref, c_ref);
%         plot(cx_ref(I_ref(c, i)), cy_ref(I_ref(c, i)), 'o')
%         contour(REF.R, REF.Z, REF.Flux,'LineColor', 'm', 'LineWidth', 2);
%     end
% end

% I_ref = reshape(I_ref, [1, 2*n]);

% zhalf = round(size(psi, 1) / 2)
% for c = 1:2
%     if c == 1
%         range = 1:zhalf;
%     else
%         range = zhalf+1:size(psi, 1);
%     end

%     for i = 1:n
%         theta = 2 * pi / n * (i - 1) + pi / n;
%         r = l * sin(theta); z = l * cos(theta);
%         line_r = [R0(c), R0(c) + r]; line_z = [Z0(c), Z0(c) + z];
%         % ある直線が作るpsiの断面を取得
%         [cx,cy,cc] = improfile(CCR, CCZ(range), psi(range, :), line_r, line_z, s, 'bilinear');
%         cc(1:ignore) = 0.01;
%         [m, J] = min(abs(cc));
%         I(c, i) = J;
%         if CONFIG.ShowFig
%             plot([line_r], [line_z], 'r');
%         end
%         plot(cx(I(c, i)), cy(I(c, i)), 'o')
%     end
% end

% v = linspace(-20, 20, 101);
% contour(CCR, CCZ, psi*1000, v)
% I= reshape(I, [1, 2*n]);