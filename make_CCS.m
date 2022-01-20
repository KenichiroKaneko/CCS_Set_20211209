function CCSDAT = make_CCS(PARAM)

    if PARAM.CCS == 0
        CCSDAT.NCCS = [];
        CCSDAT.NCCN = [];
        CCSDAT.RCCS = [];
        CCSDAT.ZCCS = [];
        CCSDAT.RCCN = [];
        CCSDAT.ZCCN = [];
        CCSDAT.RGI = [];
        CCSDAT.ZGI = [];
        return
    end
    PARAM.NE = PARAM.MSEC(1:PARAM.CCS, 1) + PARAM.MSEC(1:PARAM.CCS, 2) + PARAM.MSEC(1:PARAM.CCS, 3);

    if (PARAM.IDECCS == 0)
        % 楕円型のCCS
        CCSDAT.NCCS = PARAM.NE * 2;
        CCSDAT.NCCN = PARAM.NE * 3;

        for i = 1:PARAM.CCS
            DTHETA = 2.0 * pi / CCSDAT.NCCS(i);
            TRIG = 0.0; %三角度

            j = 1:CCSDAT.NCCS(i);
            THETA = pi / 2.0 - DTHETA * (j - 1); % CCSでの積分は、時計回り
            CCSDAT.RCCS(i, j) = PARAM.R0(i) + PARAM.RR(i) * cos(THETA + asin(TRIG) * sin(THETA));
            CCSDAT.ZCCS(i, j) = PARAM.Z0(i) + PARAM.CAPPER(i) * PARAM.RR(i) * sin(THETA);

        end

    else
        % D型のCCS
        SR(1:PARAM.CCS, 1) = PARAM.R0 - PARAM.RR;
        SR(1:PARAM.CCS, 3) = SR(1:PARAM.CCS, 1);
        SR(1:PARAM.CCS, 2) = PARAM.R0 + PARAM.RR;
        SZ(1:PARAM.CCS, 1) = PARAM.Z0 + PARAM.RR .* PARAM.CAPPER;
        SZ(1:PARAM.CCS, 2) = PARAM.Z0;
        SZ(1:PARAM.CCS, 3) = PARAM.Z0 - PARAM.RR .* PARAM.CAPPER;

        for i = 1:PARAM.CCS
            %% 曲線部のメッシュ点を定義
            II = 0;

            for j = 1:2
                M2 = PARAM.MSEC(i, j) * 2; % メッシュ点の数
                DEL = 1.0 / M2;
                GISTAT = -1.0 + (j - 1); % グザイ -1 0 1

                for k = 1:M2
                    GI = GISTAT + DEL * (k - 1); %グザイ を境界要素分に分割している
                    % Compute the values of the shape functions at the integration points
                    F1 = GI * (GI - 1.0) / 2;
                    F2 = 1.0 - GI^2;
                    F3 = GI * (GI + 1.0) / 2;
                    II = II + 1;
                    CCSDAT.RCCS(i, II) = SR(i, 1) * F1 + SR(i, 2) * F2 + SR(i, 3) * F3;
                    CCSDAT.ZCCS(i, II) = SZ(i, 1) * F1 + SZ(i, 2) * F2 + SZ(i, 3) * F3;
                end

            end

            %% 直線部のメッシュ点を定義
            II = II + 1;
            CCSDAT.RCCS(i, II) = SR(i, 3);
            CCSDAT.ZCCS(i, II) = SZ(i, 3);
            KADO = II;
            M2 = PARAM.MSEC(i, 3) * 2;
            DELR = (SR(i, 1) - SR(i, 3)) / M2;
            DELZ = (SZ(i, 1) - SZ(i, 3)) / M2;

            for k = 1:M2
                II = II + 1;
                CCSDAT.RCCS(i, II) = CCSDAT.RCCS(i, II - 1) + DELR;
                CCSDAT.ZCCS(i, II) = CCSDAT.ZCCS(i, II - 1) + DELZ;
            end

            CCSDAT.NCCS(i) = II - 1;
            CCSDAT.NCCN(i) = PARAM.NE(i) * 3;

            %% CCS上の非適合要素節点座標の作成
            CCSDAT.RCCS(i, CCSDAT.NCCS + 1) = CCSDAT.RCCS(i, 1);
            CCSDAT.ZCCS(i, CCSDAT.NCCS + 1) = CCSDAT.ZCCS(i, 1);
            CCSDAT.NCCN(i) = PARAM.NE(i) * 3; % 非適合要素
            j = 1:PARAM.NE(i);
            CCSDAT.RCCN(i, 3 * j - 2) = (5 * CCSDAT.RCCS(i, 2 * j - 1) + 5 .* CCSDAT.RCCS(i, 2 * j) - CCSDAT.RCCS(i, 2 * j + 1)) / 9;
            CCSDAT.ZCCN(i, 3 * j - 2) = (5 * CCSDAT.ZCCS(i, 2 * j - 1) + 5 .* CCSDAT.ZCCS(i, 2 * j) - CCSDAT.ZCCS(i, 2 * j + 1)) / 9;
            CCSDAT.RCCN(i, 3 * j - 1) = CCSDAT.RCCS(i, 2 * j);
            CCSDAT.ZCCN(i, 3 * j - 1) = CCSDAT.ZCCS(i, 2 * j);
            CCSDAT.RCCN(i, 3 * j) = (5 * CCSDAT.RCCS(i, 2 + j + 1) + 5 .* CCSDAT.RCCS(i, 2 * j) - CCSDAT.RCCS(i, 2 * j - 1)) / 9;
            CCSDAT.ZCCN(i, 3 * j) = (5 * CCSDAT.ZCCS(i, 2 * j + 1) + 5 .* CCSDAT.ZCCS(i, 2 * j) - CCSDAT.ZCCS(i, 2 * j - 1)) / 9;
        end

    end

    count = 1;

    for i = 1:PARAM.CCS
        CCSDAT.RCCS(i, CCSDAT.NCCS + 1) = CCSDAT.RCCS(i, 1);
        CCSDAT.ZCCS(i, CCSDAT.NCCS + 1) = CCSDAT.ZCCS(i, 1);

        I = 1:PARAM.NE(i);
        CCSDAT.RCCN(i, 3 * I - 2) = (5 .* CCSDAT.RCCS(i, 2 * I - 1) + 5 .* CCSDAT.RCCS(i, 2 * I) - CCSDAT.RCCS(i, 2 * I + 1)) / 9;
        CCSDAT.ZCCN(i, 3 * I - 2) = (5 .* CCSDAT.ZCCS(i, 2 * I - 1) + 5 .* CCSDAT.ZCCS(i, 2 * I) - CCSDAT.ZCCS(i, 2 * I + 1)) / 9;
        CCSDAT.RCCN(i, 3 * I - 1) = CCSDAT.RCCS(i, 2 * I);
        CCSDAT.ZCCN(i, 3 * I - 1) = CCSDAT.ZCCS(i, 2 * I);
        CCSDAT.RCCN(i, 3 * I) = (5 .* CCSDAT.RCCS(i, 2 * I + 1) + 5 .* CCSDAT.RCCS(i, 2 * I) - CCSDAT.RCCS(i, 2 * I - 1)) / 9;
        CCSDAT.ZCCN(i, 3 * I) = (5 .* CCSDAT.ZCCS(i, 2 * I + 1) + 5 .* CCSDAT.ZCCS(i, 2 * I) - CCSDAT.ZCCS(i, 2 * I - 1)) / 9;
    end
    

    %  Draw a fine curve to express the shape of CCS
    MAXG = 200;
    GMAX = MAXG;
    DEL = 2.0 / GMAX;

    for I = 1:PARAM.NE(i)
        for J = 1:MAXG + 1
            CJM1 = J - 1;
            GII = -1.0 + DEL * CJM1;
            F1 = GII * (GII - 1.0D0) * 0.5;
            F2 = 1.0D0 - GII^2;
            F3 = GII * (GII + 1.0D0) * 0.5;
            RGI(count)  = CCSDAT.RCCS(i, 2 * I - 1) * F1 + ...
            CCSDAT.RCCS(i, 2 * I) * F2 + CCSDAT.RCCS(i, 2 * I + 1) * F3;
            ZGI(count) = CCSDAT.ZCCS(i, 2 * I - 1) * F1 + ...
            CCSDAT.ZCCS(i, 2 * I) * F2 + CCSDAT.ZCCS(i, 2 * I + 1) * F3;
            count = count + 1;
        end

    end

    CCSDAT.RGI = RGI;
    CCSDAT.ZGI = ZGI;

end
