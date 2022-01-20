% 新たなLCFSの評価
clear all;
close all;
load('vars_exp')

figure()
hold on
% 磁気軸
[m r0] = min(min(REF.Flux));
[m z0] = min(REF.Flux(:, r0));
r0 = REF.R(r0);
z0 = REF.Z(z0);
plot(r0, z0, 'o')
contour(CCR, CCZ, psi)
contour(REF.R, REF.Z ,REF.Flux)

l = 0.6; % 直線の長さ
s = 100; % 直線の内部を補完する数
ignore = round(s * 0.3) % 内側から無視する数(中心付近は磁束が乱れる)
n = 32; % 放射状にいくつの線を引くか

% ReferenceのLCFSの評価
for i = 1:n
    theta = 2 * pi / n * (i - 1) + pi / n;
    r = l * sin(theta);
    z = l * cos(theta);
    line_r = [r0, r0 + r];
    line_z = [z0, z0 + z];
    % ある直線が作るpsiの断面を取得
    [cx_ref, cy_ref, c_ref] = improfile(REF.R, REF.Z, REF.Flux, line_r, line_z, s, 'bilinear');
    c_ref(1:ignore) = [];
    % LCFS(psi=0)との交点を探索
    [m, J] = min(abs(c_ref));
    I_ref(i) = J + ignore;
    plot3(cx_ref, cy_ref, c_ref);
    plot(cx_ref(I_ref(i)), cy_ref(I_ref(i)), 'o')
    plot(cx_ref(ignore), cy_ref(ignore), 'o')
    
end

for i = 1:n
    theta = 2 * pi / n * (i - 1) + pi / n;
    r = l * sin(theta);
    z = l * cos(theta);
    line_r = [r0, r0 + r];
    line_z = [z0, z0 + z];
    % ある直線が作るpsiの断面を取得
    [cx,cy,c] = improfile(CCR, CCZ, psi, line_r, line_z, s, 'bilinear');
    c(1:ignore) = [];
    % LCFS(psi=0)との交点を探索
    [m, J] = min(abs(c));
    I(i) = J + ignore;
    plot3(cx,cy,c);
    plot(cx(I(i)), cy(I(i)), 'o')
    plot(cx(ignore), cy(ignore), 'o')
end

% matlabは1 origin なのでインデックスを-1する
% 正しい距離を計算するのも、インデックスで計算するのも結果は同じ
1 / n * sum(abs(I_ref - I) ./ (I_ref-1)) * 100