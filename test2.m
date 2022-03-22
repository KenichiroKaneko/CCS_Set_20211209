clear all;
close all;

% f1 = figure()

load("cal/result_2033_良好な1段目との比較/errs")

for i = 1:25

    disp(i)
    % load("cal/result_2033_良好なKUPを探索誤差なし" + num2str(i))
    % load("cal/result_2033_良好な1段目との比較/"+num2str(i))
    load("cal/result_2033_1段目でLを検証/"+num2str(i))
    REF1 = load('vars_result_良い', "CCR", "CCZ");
    load('psi_add', "psi_add", "CCR", "CCZ"); 
    psiBig = imresize(psi, size(REF.Flux));
    psi = psiBig + psi_add;
    CCR = REF.R;
    CCZ = REF.Z;
    % v = linspace(-20, 20, 101);

    % figure(f1)
    % subplot(1,7,i-18)
    % contour(CCR, CCZ, psi*1000, v)

    % figure()
    % v = linspace(-20, 20, 101);
    % hold on
    % psiBig = psiBig + psi_add;
    % contour(REF.R, REF.Z, psiBig*1000, v, 'k');
    % contour(REF.R, REF.Z, psiBig, [0 0], 'c', 'LineWidth', 2)
    % contour(REF.R, REF.Z, REF.Flux*1000, v, 'm--'); % 正解;
    % contour (REF.R, REF.Z, REF.Flux, [0 0], 'm', 'LineWidth', 2);
    % plot(CCSDAT.RGI, CCSDAT.ZGI, CCSDAT.RGI, -CCSDAT.ZGI, 'color','y','LineWidth', 2);
    % plot(WALL.RWALL_list, WALL.ZWALL_list, '-k'); % 容器壁 VacuumVesselMeshPoints
    % hold off
    % xlabel({'r (m)'});
    % ylabel({'z (m)'});
    % title("Reconstructed flux")
    % axis equal

    err_LCFS(i) = evaluate_LCFS(psi, REF, PARAM, CONFIG, CCR, CCZ, 0);
    % FF-AA*X 残差
    residual(i) = norm(AA * FFout - FF);

end


[m, I] = min(residual)
[m, II] = min(err_LCFS)

figure()
hold on
xlabel("LCFS誤差")
ylabel("残差 ||Dp^* - g||")
plot(err_LCFS, residual, 'o')
plot(err_LCFS(I), residual(I), 'c*')
plot(err_LCFS(II), residual(II), 'm*')


for i = 1:length(err_LCFS)
    text(err_LCFS(i), (residual(i)), ['\leftarrow', num2str(i + 49)])
end

figure()
hold on
plot(residual)
plot(I, residual(I), 'o')

figure()
plot(err_LCFS)
