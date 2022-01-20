function cal_MRE(PARAM, CONFIG, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, A, ss, vv, uu, X, FC, FF)
    % 各センサーの相対誤差を計算
    % センサーの順番はprobe,Flux,CCSだが、コードでは間違えている
    % 今度直す

    NFLX = SENSOR_FLXLP.NUM;
    NAPB = SENSOR_TPRB.NUM + SENSOR_NPRB.NUM;
    NCCN = sum(CCSDAT.NCCN);

    NFLXl = 1;
    NFLXr = NFLX;
    NAPBl = 1 + NFLXr;
    NAPBr = NFLXr + NAPB;
    NCCSl = 1 + NAPBr;
    NCCSr = NAPBr + NCCN;

    indexB_l = 1;
    indexB_r = NAPB;
    indexFL_l = 1 + indexB_r;
    indexFL_r = indexFL_l + NFLX - 1;
    indexCCS_l = 1 + indexFL_r;
    indexCCS_r = indexCCS_l + NCCN - 1;

    % RE = ((A * X - FC') ./ FC').^2;
    RE = abs(A * X - FF');

    MREAPB = norm(RE(indexB_l:indexB_r)) / NAPB;
    MREFlux = norm(RE(indexFL_l:indexFL_r)) / NFLX;
    MRECCS = norm(RE(indexCCS_l:indexCCS_r)) / NCCN;

    MREs = [MREAPB, MREFlux, MRECCS];
    % disp(MREs);

    if CONFIG.ShowFig
        txt1 = sprintf("Probe Sensor  Num:%4d  MRE:%.3e\n", NAPB, MREAPB);
        txt2 = sprintf("Flux  Sensor  Num:%4d  MRE:%.3e\n", NFLX, MREFlux);
        txt3 = sprintf("CCS Node      Num:%4d  MRE:%.3e\n", NCCN, MRECCS);
        disp(txt1 + txt2 + txt3)
    end


    % % ふし目を探索
    % R = SENSOR_FLXLP.R;
    % Z = SENSOR_FLXLP.Z;
    % count_FLUX = 0;

    % for i = 1:SENSOR_FLXLP.NUM - 1

    %     if (R(i) ~= R(i + 1) && Z(i) ~= Z(i + 1)) || (i == 1) || (i == NAPB)
    %         count_FLUX = count_FLUX + 1;
    %         corner_FLUX(count_FLUX) = i;
    %     end

    % end

    % R = SENSOR_NPRB.R;
    % Z = SENSOR_NPRB.Z;
    % count_NPRB = 0;

    % for i = 1:SENSOR_NPRB.NUM - 1

    %     if (R(i) ~= R(i + 1) && Z(i) ~= Z(i + 1)) || (i == 1) || (i == SENSOR_NPRB.NUM - 1)
    %         count_NPRB = count_NPRB + 1;
    %         corner_NPRB(count_NPRB) = NFLX + i
    %     end

    % end

    % plot(NFLXl:NFLXr, FF(NFLXl:NFLXr))

    if CONFIG.ShowFig
        AX = A * X;
        figure(); title('各センサー位置の再構成結果');

        subplot(3, 1, 1);
        hold on
        plot(1:NAPB, AX(indexB_l:indexB_r), '-o');
        plot(1:NAPB, FF(indexB_l:indexB_r), '-o');
        title('sensor B'); legend('A*X', 'FF', 'location', 'eastoutside');
        yline(0);
        hold off

        subplot(3, 1, 2);
        hold on
        plot(1:NFLX, AX(indexFL_l:indexFL_r), '-o');
        plot(1:NFLX, FF(indexFL_l:indexFL_r), '-o');
        title('sensor Flux Loop'); legend('A*X', 'FF', 'location', 'eastoutside');
        yline(0);
        hold off

        subplot(3, 1, 3);
        hold on
        plot(1:NCCN, AX(indexCCS_l:indexCCS_r), '-o');
        plot(1:NCCN, FF(indexCCS_l:indexCCS_r), '-o');
        title('CCS node'); legend('A*X', 'FF', 'location', 'eastoutside');
        hold off

        % x = [NFLXr, NFLXr, NAPBl, NAPBr, NCCSl, NCCSr];
        % scatter(x, FF(x), "m*")

        % for i = 1:count_FLUX
        %     scatter(corner_FLUX(i), FC(corner_FLUX(i)), 'm*')
        % end

        % for i = 1:count_NPRB
        %     scatter(corner_NPRB(i), FC(corner_NPRB(i)), 'm*')
        % end

        figure()
        plot(NFLXl:NCCSr, A * X - FF', '-o');
        % hold on
        % x = [NFLXr, NFLXr, NAPBl, NAPBr, NCCSl, NCCSr];
        % scatter(x, zeros(1, length(x)), "m*")
        title("A * X'-FF");
        % legend();
    end

    % save("vars_inMRE  ");
    % error('f', A1)
end
