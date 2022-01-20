function [PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, POS] = make_sensor(PARAM, CONFIG, BZ, BR, FL, POS)
    %% 再構成する時間でのデータと位置や角度をchごとに対応させる
    senpos_B = fileread(PARAM.temporary_file_directory + PARAM.MP_pos);
    senpos_B = strsplit(senpos_B, {'\n', '\t', '\r'});
    senpos_FL = fileread(PARAM.temporary_file_directory + PARAM.FL_pos);
    senpos_FL = strsplit(senpos_FL, {'\n', '\t', '\r'});
    SENSOR_TPRB.NUM = floor((length(senpos_B) - 5) / 5);
    SENSOR_FLXLP.NUM = floor((length(senpos_FL) - 5) / 5);
    if CONFIG.DataType == 'exp'
        SENSOR_TPRB.R = str2double(senpos_B(6:5:5 * SENSOR_TPRB.NUM + 5));
        SENSOR_TPRB.Z = str2double(senpos_B(7:5:5 * SENSOR_TPRB.NUM + 5));
        SENSOR_FLXLP.R = str2double(senpos_FL(6:5:5 * SENSOR_FLXLP.NUM + 5));
        SENSOR_FLXLP.Z = str2double(senpos_FL(7:5:5 * SENSOR_FLXLP.NUM + 5));
    elseif CONFIG.DataType == 'sol'
        SENSOR_TPRB.R = POS.Br;
        SENSOR_TPRB.Z = POS.Bz;
        SENSOR_FLXLP.R = POS.FLr;
        SENSOR_FLXLP.Z = POS.FLz;
        SENSOR_FLXLP.NUM = length(POS.FLz);
        SENSOR_TPRB.NUM = length(POS.Br);
    end

    % EFがつくる磁場、磁束を差し引く
    EF_voltage = 120;
    [Bz_EF_at_sensor_f, Psi_EF_at_sensor_f] = EF_calc_for_CCS_probe(SENSOR_FLXLP.R, SENSOR_FLXLP.Z, 1, 35, EF_voltage);
    [Bz_EF_at_sensor_b, Psi_EF_at_sensor_b] = EF_calc_for_CCS_probe(SENSOR_TPRB.R, SENSOR_TPRB.Z, 1, 40, EF_voltage);
    
    % BZ
    SENSOR_TPRB.TET = atan2(BZ, BR); % Theta センサーの角度
    
    if CONFIG.DataType == 'exp'
        BZ = BZ + Bz_EF_at_sensor_b;
        for i = 1:SENSOR_TPRB.NUM
            if SENSOR_TPRB.R(i) > 0.5
                SENSOR_TPRB.TET(i) = -pi / 2;
            else
                SENSOR_TPRB.TET(i) = pi / 2;
            end
        end
    end
    
    SENSOR_TPRB.TPRB = sqrt(BR.^2 + BZ.^2);
    SENSOR_TPRB.ITYPE = 1; % 使っていない

    %% No NPRB
    SENSOR_NPRB.NUM = 0;
    SENSOR_NPRB.R = [];
    SENSOR_NPRB.Z = [];
    SENSOR_NPRB.TET = [];
    SENSOR_NPRB.NPRB = [];
    SENSOR_NPRB.ITYPE = [];

    % FL
    if CONFIG.DataType == 'exp'
        FL = FL + Psi_EF_at_sensor_f;
    end

    if CONFIG.DevideFlux
        SENSOR_FLXLP.FLXLP = FL ./ (pi .* SENSOR_FLXLP.R .* SENSOR_FLXLP.R);
    else
        SENSOR_FLXLP.FLXLP = FL;
    end

    SENSOR_FLXLP.TET = ones(1, SENSOR_FLXLP.NUM);
    SENSOR_FLXLP.ITYPE = zeros(1, SENSOR_FLXLP.NUM);

    % CCSセンター位置の決定
    if CONFIG.AutoCCSZPos
        RR = SENSOR_TPRB.R; ZZ = SENSOR_TPRB.Z;
        B = SENSOR_TPRB.TPRB;
        % RR = SENSOR_FLXLP.R; ZZ = SENSOR_FLXLP.Z;
        % B = SENSOR_FLXLP.FLXLP;
        RR0index = RR < min(RR) + 0.1;
        B = B(RR0index); ZZ = ZZ(RR0index); RR = RR(RR0index);
        % ZZ0index = and(ZZ < 0.4, ZZ > -0.4);
        % RR = RR(ZZ0index); ZZ = ZZ(ZZ0index); B = B(ZZ0index);
        ZZ = ZZ'; B = B';
        f = fit(ZZ, B, 'smoothingspline', 'SmoothingParam', 1);
        x = fnzeros(fnder(f.p)); % 微分して勾配が０の点を探している
        x = unique(x(:));
        j = 1;
        lo_max = 0;

        for i = 1:length(x)
            % 磁気軸のありそうな範囲の中で極大値を探索
            if (fnval(fnder(f.p), x(i) - 1e-5) > 0 & fnval(fnder(f.p), x(i) + 1e-5) < 0 & abs(x(i)) < 0.45)
                lo_max(j) = x(i);
                j = j + 1;
            end

        end

        % x = sort(lo_max)
            x = sort(lo_max, 'descend', 'ComparisonMethod', 'abs');
        lmax = islocalmax(B);

        if CONFIG.ShowFig
            figure()
            hold on
            plot(ZZ, B);
            fnplt(f.p);
            plot(x, fnval(f.p, x), 'o');
            xlim([-1 1]);
            plot(ZZ(lmax), B(lmax), 'b*');
            view(90, 90);
            ax = gca;
            ax.XDir = 'reverse';
        end

        if length(x) > 1

            if (abs(x(1)) + abs(x(2))) > 0.25
                PARAM.CCS = 2;
                kanekoZ = [0.22, -0.28];
                % kanekoZ = [0.06, -0.06];

                for i = 1:2
                    % PARAM.R0(i) = 0.35;
                    PARAM.Z0(i) = x(i);
                    % PARAM.Z0(i) = kanekoZ(i);
                    PARAM.RR(i) = PARAM.RR(1);
                    PARAM.CAPPER(i) = PARAM.CAPPER(1);
                    PARAM.NE(i) = PARAM.NE(1);
                    PARAM.MSEC(i, 1) = PARAM.MSEC(1, 1);
                    PARAM.MSEC(i, 2) = PARAM.MSEC(1, 2);
                    PARAM.MSEC(i, 3) = PARAM.MSEC(1, 3);
                    PARAM.SOU(i) = PARAM.SOU(1);

                end

            end

        end

        % CCS面の自動決定R方向
        if CONFIG.AutoCCSRPos
            PARAM.R0 = CalcPlasmaCenter(PARAM, CONFIG, PARAM.Z0);
        end

    end

    % 死んでそうな磁場センサーを消す
    i = PARAM.dead_BZ;
    SENSOR_TPRB.R(i) = [];
    SENSOR_TPRB.Z(i) = [];
    SENSOR_TPRB.TET(i) = [];
    SENSOR_TPRB.TPRB(i) = [];
    SENSOR_TPRB.NUM = length(SENSOR_TPRB.TPRB);
    if CONFIG.DataType == 'sol'
        POS.Br(i) = [];
        POS.Bz(i) = [];
    end
    % 死んでそうなFLXLPセンサーを消す
    i = PARAM.dead_FL;
    SENSOR_FLXLP.R(i) = [];
    SENSOR_FLXLP.Z(i) = [];
    SENSOR_FLXLP.TET(i) = [];
    SENSOR_FLXLP.FLXLP(i) = [];
    SENSOR_FLXLP.NUM = length(SENSOR_FLXLP.FLXLP);
    if CONFIG.DataType == 'sol'
        POS.FLr(i) = [];
        POS.FLz(i) = [];
    end
end


function [Bz, Psi] = EF_calc_for_CCS_probe(r, z, r_size, z_size, V)
    R_EF = 855 * 1e-3;
    z_EF = 1.05;
    Turn_EF = 200;
    mu0 = 4 * pi * 1e-7;
    V = 120;
    I =- (0.849 * (1.19 * V - 5.32) - 5.56);
    Psi = zeros(z_size, r_size);

    for i = 1:1
        alpha_u = sqrt((R_EF + r).^2 + (z - z_EF).^2);
        alpha_l = sqrt((R_EF + r).^2 + (z + z_EF).^2);
        k2_u = (4 * R_EF * r) ./ alpha_u.^2;
        k_u = sqrt(k2_u);
        k2_l = (4 * R_EF * r) ./ alpha_l.^2;
        k_l = sqrt(k2_l);
        m_u = k2_u;
        m_l = k2_l;
        [K_u, E_u] = ellipke(m_u);
        [K_l, E_l] = ellipke(m_l);
        Bz_upper = ((mu0 * Turn_EF * I) ./ (2 * pi)) .* (1 ./ alpha_u) .* (K_u + ((R_EF^2 - r.^2 - (z - z_EF).^2) ./ ((R_EF - r).^2 + (z - z_EF).^2)) .* E_u);
        Bz_lower = ((mu0 * Turn_EF * I) ./ (2 * pi)) .* (1 ./ alpha_l) .* (K_l + ((R_EF^2 - r.^2 - (z + z_EF).^2) ./ ((R_EF - r).^2 + (z + z_EF).^2)) .* E_l);
        Bz = Bz_upper + Bz_lower;

        A_phai_u = (mu0 * Turn_EF * I) ./ (pi) .* sqrt(R_EF ./ r) .* (1 ./ k_u) .* ((1 - k_u.^2 ./ 2) .* K_u - E_u);
        A_phai_l = (mu0 * Turn_EF * I) ./ (pi) .* sqrt(R_EF ./ r) .* (1 ./ k_l) .* ((1 - k_l.^2 ./ 2) .* K_l - E_l);

        psi_u = 2 * pi * r .* A_phai_u;
        psi_l = 2 * pi * r .* A_phai_l;
        Psi = psi_u + psi_l;

    end

end
