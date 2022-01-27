function FF = make_FF(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT)

    % figure() 
    % PARAM.SIGM = 0.03;
    % hold on
    % plot([SENSOR_TPRB.TPRB, SENSOR_FLXLP.FLXLP])

    rng(PARAM.SEED)
    if max(SENSOR_TPRB.TPRB) > 0
        e1 = SENSOR_TPRB.TPRB .* PARAM.SIGM .* randn(1, SENSOR_TPRB.NUM);
    else
        e1 = (-1) * SENSOR_TPRB.TPRB .* PARAM.SIGM .* randn(1, SENSOR_TPRB.NUM);
    end

    if max(SENSOR_FLXLP.FLXLP) > 0
        e2 = SENSOR_FLXLP.FLXLP .* PARAM.SIGM .* randn(1, SENSOR_FLXLP.NUM);
    else
        e2 = (-1) * SENSOR_FLXLP.FLXLP .* PARAM.SIGM .* randn(1, SENSOR_FLXLP.NUM);
    end

    SENSOR_TPRB.TPRB = SENSOR_TPRB.TPRB + e1;
    SENSOR_FLXLP.FLXLP = SENSOR_FLXLP.FLXLP + e2;


    % plot([SENSOR_TPRB.TPRB, SENSOR_FLXLP.FLXLP])

    % error('error description', A1)
    FF = [SENSOR_TPRB.TPRB, SENSOR_NPRB.NPRB, (SENSOR_FLXLP.FLXLP ./ 2 ./ pi)];
    FF(end + 1:end + sum(CCSDAT.NCCN)) = 0.0;

    if 0
        %% FFDAT (w. noise)
        rng(PARAM.SEED);
        GASDEV = randn(SENSOR_TPRB.NUM + SENSOR_NPRB.NUM + SENSOR_FLXLP.NUM, 1);

        i = 1:SENSOR_TPRB.NUM + SENSOR_NPRB.NUM + SENSOR_FLXLP.NUM;
        FF(i) = FF(i) * (1.0 + PARAM.SIGM * GASDEV(i, 1));
        % 磁場センサーには磁場センサーの最大値を基準にしたノイズを載せる
        % 磁束も同様
        FF(i) = FF(i) + max(SENSOR_TPRB.TPRB);
    end

end
