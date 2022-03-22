function [FF, FC, AA, factors] = weight(CONFIG, FF, FC, AA, SENSOR_NPRB, SENSOR_TPRB, SENSOR_FLXLP, CCSDAT)
    %% INOMOTO start
    % 規格化に使うそれぞれのセンサーのインデックス
    % センサーを使わない場合でも正しく処理を行えるように小細工してある
    if SENSOR_NPRB.NUM + SENSOR_TPRB.NUM ~= 0
        ind_b = 1:SENSOR_NPRB.NUM + SENSOR_TPRB.NUM;
    else
        ind_b = 1:1;
    end

    if SENSOR_FLXLP.NUM ~= 0
        ind_FLXLP = SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + 1:...
        SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM;
    else
        ind_FLXLP = 1:1;
    end

    if CCSDAT.NCCN ~= 0
        ind_CCS = SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM + 1:...
        SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM + sum(CCSDAT.NCCN);
    else
        ind_CCS = 1:1;
    end

    % 重要なセンサーに重み付け
    % load('vars_roop_合体後1', 'err_LCFS')
    % FF(ind_b) = FF(ind_b) .* (err_LCFS(ind_b) + 1)*10;
    % AA(ind_b, :) = AA(ind_b, :) .* (err_LCFS(ind_b) + 1)'*10;
    % FF(ind_FLXLP) = FF(ind_FLXLP) .* (err_LCFS(ind_FLXLP) + 1)*10;
    % AA(ind_FLXLP, :) = AA(ind_FLXLP, :) .* (err_LCFS(ind_FLXLP) + 1)'*10;
    
    
    %% Modified_MI 20211103
    ind_B_IN = 1:18;
    ind_B_OUT = 20:39;
    ind_FLXLP_IN = [1:19] + SENSOR_TPRB.NUM + SENSOR_NPRB.NUM;
    ind_FLXLP_OUT = [21:34] + SENSOR_TPRB.NUM + SENSOR_NPRB.NUM;
    size(FF);
    %% Modified_MI 20211103 kokomade    

    switch CONFIG.Nmrz
        case 'ave'
            % 平均が１になるように規格化→センサー信号が全て０付近の場合悪くなりそう
            bfactor = 1 / mean(FF(ind_b));
            fluxfactor = 1 / mean(FF(ind_FLXLP));
            % ccsfactor = 1 / mean(FF(ind_CCS));
            ccsfactor = 1;
        case 'max'
            % 絶対値の最大値が１になるように規格化
            bfactor = 1 / max(abs(FF(ind_b)));
            fluxfactor = 1 / max(abs(FF(ind_FLXLP)));
            % ccsfactor = 1 / max(abs(FF(ind_CCS)));
            ccsfactor = 1;
        case 'norm'
            % ノルムで割る規格化
            bfactor = 1 / norm(FF(ind_b), 1);
            fluxfactor = 1 / norm(FF(ind_FLXLP), 1);
            ccsfactor = 1;
            % ccsfactor = 1 / norm(FF(ind_CCS), 1);
        case '10'
            % 手動の規格化
            bfactor = 1;
            fluxfactor = 10;
            ccsfactor = 1;
        case '100'
            % 手動の規格化
            bfactor = 1;
            fluxfactor = 100;
            ccsfactor = 1;
        case 'no'
            % 規格化なし
            bfactor = 1;
            fluxfactor = 1;
            ccsfactor = 1;
        otherwise    
            % 規格化なし
            bfactor = 1;
            fluxfactor = 1;
            ccsfactor = 1;
    end

    % % fluxのノルムをbと同じにする規格化
    % bfactor = 1;
    % fluxfactor = norm(FF(ind_b), 1) / norm(FF(ind_FLXLP), 1);
    % ccsfactor = norm(FF(ind_b), 1) / norm(FF(ind_CCS), 1);

    
    if SENSOR_NPRB.NUM + SENSOR_TPRB.NUM == 0
        bfactor = 1;
    end

    if SENSOR_FLXLP.NUM == 0
        fluxfactor = 0;
    end

    if CCSDAT.NCCN == 0
        ccsfactor = 0;
    end

    %
    FF(ind_b) = FF(ind_b) * bfactor;
    AA(ind_b, :) = AA(ind_b, :) * bfactor;
    % FC(ind_b) = FC(ind_b) * bfactor;
    %% Modified_MI 20211103
    FF(ind_FLXLP) = FF(ind_FLXLP) * fluxfactor;
    AA(ind_FLXLP, :) = AA(ind_FLXLP, :) * fluxfactor;
    % FC(ind_FLXLP) = FC(ind_FLXLP) * fluxfactor;
    % FF(ind_FLXLP_IN) = FF(ind_FLXLP_IN) * fluxfactor;
    % AA(ind_FLXLP_IN, :) = AA(ind_FLXLP_IN, :) * fluxfactor;
    % FC(ind_FLXLP_IN) = FC(ind_FLXLP_IN) * fluxfactor;
    % FLXLP = [FFDAT(ind_FLXLP_IN) * fluxfactor FFDAT(ind_FLXLP_OUT)];
    %% Modified_MI 20211103 kokomade
    FF(ind_CCS) = FF(ind_CCS) * ccsfactor;
    AA(ind_CCS, :) = AA(ind_CCS, :) * ccsfactor;
    % FC(ind_CCS) = FC(ind_CCS) * ccsfactor;
    factors = [bfactor, fluxfactor, ccsfactor];



end