function show_results(PARAM, CONFIG, CCR, CCZ, CCSDAT, REF, psi, SENSOR_TPRB, SENSOR_FLXLP, WALL,...
    ExtCOIL, DISF, t, Ip, FFout)

    % FFoutのプロット
    figure('Name', 'FFout')
    hold on
    plot(FFout, 'o-')
    n = sum(CCSDAT.NCCN) * 2;
    xline(n, 'r--')
    title('CCS / eddy current')
    
    % コイル電流などの表示
    if CONFIG.DataType == 'exp'
        % xlabel('Time [ms]');
        % ylabel('Current [kA]');
        figure('Name', 'Plasma/Coil Current')
        time = (1:30000)*5e-4;
        subplot(5,1,1)
        hold on
        plot(time,Ip,'b','DisplayName','I plasma')
        plot(t,Ip(t),'o','DisplayName','Iccs plasma')
        xlim([8.5, 10.5])
        ylim([-100, 150])
        legend({}, 'Location', 'eastoutside','FontSize',14)
        for i = 1:4
            j = i*2-1;
            subplot(5,1,i+1)
            hold on
            plot(time, ExtCOIL.I_sig(j,:), 'r', 'DisplayName', ['Smoothed PF #' num2str(i) 'L'])
            plot(time, ExtCOIL.I_sig(j+1,:), 'b', 'DisplayName', ['Smoothed PF #' num2str(i) 'U'])
            xlim([8.5, 10.5])
            ylim([-60, 70])
            legend({}, 'Location', 'eastoutside','FontSize',14)
        end
    end

    % 渦電流再構成結果の表示
    figure('Name', 'Eddy Current Plofile', 'NumberTitle', 'off')
    hold on
    plot(DISF.dis, DISF.jeddy, '-ko', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 2)

    for i = 1:length(WALL.WALL_dist_list)
        xline(WALL.WALL_dist_list(i), '--r');
        % xline(WALL.WALL_dist_list(i)-WALL.WALL_dist_list(2), '--r');
    end

    xlabel({'Distance (m)'});
    ylabel({'Eddy Current Density (MA/m)'});

    % センサー位置と番号の表示
    figure('Name', 'Sensor position')
    subplot(1,2,1)
    hold on
    % t = 1:SENSOR_TPRB.NUM; t = num2str(t');
    % text(SENSOR_TPRB.R, SENSOR_TPRB.Z, t)
    plot(SENSOR_TPRB.R, SENSOR_TPRB.Z, 'ro')
    plot(WALL.RWALL_list, WALL.ZWALL_list, '-k'); % 容器壁 VacuumVesselMeshPoints
    title('Magnetic probe')
    axis equal
    subplot(1,2,2)
    hold on
    % t = 1:SENSOR_FLXLP.NUM; t = num2str(t');
    % text(SENSOR_FLXLP.R, SENSOR_FLXLP.Z, t)
    plot(SENSOR_FLXLP.R, SENSOR_FLXLP.Z, 'bo')
    plot(WALL.RWALL_list, WALL.ZWALL_list, '-k'); % 容器壁 VacuumVesselMeshPoints
    title('Flux loop')
    axis equal


    % 再構成結果の表示
    % psiBig = imresize(psi, size(REF.Flux));
    v = linspace(-20, 20, 101);
    figure()
    subplot(1, 2, 1);
    hold on
    % contour(REF.R, REF.Z, psiBig*1000, v); % 再構成結果
    % contour(REF.R, REF.Z, psiBig, [0 0], 'k', 'LineWidth', 2); % LCFS
    contour(CCR, CCZ, psi*1000, v); % 再構成結果
    contour(CCR, CCZ, psi, [0 0], 'k', 'LineWidth', 2); % LCFS
    plot(SENSOR_TPRB.R, SENSOR_TPRB.Z, 'ro')
    plot(SENSOR_FLXLP.R, SENSOR_FLXLP.Z, 'bo')
    plot(WALL.REV, WALL.ZEV, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8) % VacuumVesselNodePoints
    plot(WALL.RSEC, WALL.ZSEC, 'mo', 'MarkerSize', 16) % VacuumVesselSegments
    plot(WALL.RWALL_list, WALL.ZWALL_list, '-k'); % 容器壁 VacuumVesselMeshPoints
    plot(CCSDAT.RGI, CCSDAT.ZGI, CCSDAT.RGI, -CCSDAT.ZGI, 'color','y','LineWidth', 2);
    if CONFIG.DataType =='sol'
        contour(REF.R, REF.Z, REF.Flux, [0 0], 'LineColor', 'm', 'LineWidth', 2);
    end
    hold off
    xlabel({'r (m)'});
    ylabel({'z (m)'});
    title("Reconstructed flux")
    axis equal

    subplot(1, 2, 2);
    hold on
    if CONFIG.DataType =='sol'
        % REFpsi = imresize(REF.Flux, [size(psi)]);
        contour(REF.R, REF.Z, REF.Flux*1000, v, 'm'); % 正解;
        contour (REF.R, REF.Z, REF.Flux, [0 0], 'm', 'LineWidth', 2);
    elseif CONFIG.DataType == 'exp'
        % contour(mag.rr, mag.zz, squeeze(mag.data_line(PARAM.time_step, :, :)), 'LineColor', 'm')
        % contour(mag.rr, mag.zz, squeeze(mag.data_line(PARAM.time_step, :, :)), [-1e6 0 1e6], 'm', 'Linewidth', 2)
        % legend('Reconstructed Flux', 'Reconstructed LCFS');
    end
    % contour(REF.R, REF.Z, psiBig*1000, v, 'k');
    % contour(REF.R, REF.Z, psiBig, [0 0], 'c', 'LineWidth', 2)
    contour(CCR, CCZ, psi*1000, v,'k');
    contour(CCR, CCZ, psi, [0 0], 'LineColor', 'c', 'LineWidth', 2);
    plot(CCSDAT.RGI, CCSDAT.ZGI, CCSDAT.RGI, -CCSDAT.ZGI, 'color','y','LineWidth', 2);
    plot(WALL.RWALL_list, WALL.ZWALL_list, '-k'); % 容器壁 VacuumVesselMeshPoints
    hold off
    xlabel({'r (m)'});
    ylabel({'z (m)'});
    title("Reconstructed flux")
    axis equal
    % var = load('vars_CCS_temp');
    % figure()
    % hold on
    % v = linspace(-20, 20, 101);
    % contour(REF.R, REF.Z, REF.Flux*1000, v, 'k'); % 正解;
    % contour (REF.R, REF.Z, REF.Flux, [0 0], 'k', 'LineWidth', 2, 'DisplayName', 'Reference');
    % contour(CCR, CCZ, psi, [0 0], 'm', 'LineWidth', 2, 'DisplayName', '0% Noise');
    % contour(CCR, CCZ, var.psi, [0 0], 'c', 'LineWidth', 2, 'DisplayName', '3% Noise');
    % plot(WALL.RWALL_list, WALL.ZWALL_list, '-k'); % 容器壁 VacuumVesselMeshPoints
    % % legend();
    % xlabel({'r (m)'});
    % ylabel({'z (m)'});
    % axis equal

end