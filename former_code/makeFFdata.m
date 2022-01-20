%% maksFFdata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% maksFFdata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% maksFFdata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% maksFFdata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function FFDAT = makeFFdata(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP)

    %% FFDAT0 (w/o noise)
    for i = 1:SENSOR_TPRB.NUM
        FFDAT0(i) = SENSOR_TPRB.TPRB(i);
    end

    for i = 1:SENSOR_NPRB.NUM
        FFDAT0(i + SENSOR_TPRB.NUM) = SENSOR_NPRB.NPRB(i);
    end

    for i = 1:SENSOR_FLXLP.NUM
        FFDAT0(i + SENSOR_TPRB.NUM + SENSOR_NPRB.NUM) = SENSOR_FLXLP.FLXLP(i) / 2 / pi;
    end

    if 0 % startsWith(PARAM.input_file_directory, "UTST")
        %% FFDAT (w. noise)
        rng(PARAM.SEED);
        GASDEV = randn(SENSOR_TPRB.NUM + SENSOR_NPRB.NUM + SENSOR_FLXLP.NUM, 1);

        for i = 1:SENSOR_TPRB.NUM + SENSOR_NPRB.NUM + SENSOR_FLXLP.NUM
            FFDAT(i) = FFDAT0(i) * (1.0 + PARAM.SIGM * GASDEV(i, 1));
        end

    else
        FFDAT = FFDAT0;
    end

end

%% maksFFdata kokomade!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
