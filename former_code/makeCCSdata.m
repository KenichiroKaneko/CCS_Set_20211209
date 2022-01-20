%% maksCCSdata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% maksCCSdata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% maksCCSdata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% maksCCSdata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function CCSDAT = makeCCSdata(PARAM, GHR, GHZ)

    PARAM.NE = PARAM.MSEC(1:PARAM.CCS, 1) + PARAM.MSEC(1:PARAM.CCS, 2) + PARAM.MSEC(1:PARAM.CCS, 3);
    %PARAM.NCCS = NE*2;

    ICONT = 0;
    %  9999 CONTINUE ! ICONTの更新（CCSを作り直す）!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (ICONT == 0) % ICONT>0 ならCCSを作り直す。

        if (PARAM.IDECCS > 0) %! 楕円型CCS or D型CCS
            %*******************************************************
            %  'D'型CCS                                         ****
            %*******************************************************
            CCSDAT = D_CCS(PARAM);
        else
            % *******************************************************
            %   楕円型CCS                                        ****
            % *******************************************************
            CCSDAT.NCCS = PARAM.NE * 2;
            CCSDAT.NCCN = PARAM.NE * 3;

            for i = 1:PARAM.CCS
                DTHETA = 2.0 * pi / CCSDAT.NCCS(i);
                TRIG = 0.0; %三角度

                for j = 1:CCSDAT.NCCS(i)
                    THETA = pi / 2.0 - DTHETA * (j - 1); % CCSでの積分は、時計回り
                    CCSDAT.RCCS(i, j) = PARAM.R0(i) + PARAM.RR(i) * cos(THETA + asin(TRIG) * sin(THETA));
                    CCSDAT.ZCCS(i, j) = PARAM.Z0(i) + PARAM.CAPPER(i) * PARAM.RR(i) * sin(THETA);
                    % CCSDAT.RCCS(i, j) = PARAM.R0(i) + PARAM.CAPPER(i) * PARAM.RR(i) * cos(THETA + asin(TRIG) * sin(THETA));
                    % CCSDAT.ZCCS(i, j) = PARAM.Z0(i) + PARAM.RR(i) * sin(THETA);
                end

            end

        end

    else
        % *******************************************************
        %    LCMSに相似なCCS    (ICONT>0 のとき)             ****
        % *******************************************************
        for i = 1:PARAM.CCS

            for j = 1:CCSDAT.NCCS(i)
                k = CCSDAT.NCCS(i) + 2 - j;

                if (k > NCCS(i))
                    k = 1;
                end

                CCSDAT.RCCS(i, j) = GHR(k);
                CCSDAT.ZCCS(i, j) = GHZ(k);
            end

        end

    end

    % *******************************************************
    %  CCS上の非適合要素節点座標の作成
    % *******************************************************
    fid12 = fopen([PARAM.temporary_file_directory '/MeshPoints.txt'], 'w');
    fid13 = fopen([PARAM.temporary_file_directory '/DiscontinuousNodePoints.txt'], 'w');
    fid14 = fopen([PARAM.temporary_file_directory '/CurvCCS.txt'], 'w');
    fid15 = fopen([PARAM.temporary_file_directory '/CurvCCS_Final.txt'], 'w');

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

        for j = 1:CCSDAT.NCCS(i)
            fprintf(fid12, '%d %d\n', CCSDAT.RCCS(i, j), CCSDAT.ZCCS(i, j));
        end

        fprintf(fid12, '%d %d\n', CCSDAT.RCCS(i, 1), CCSDAT.ZCCS(i, 1));

        for j = 1:CCSDAT.NCCN(i)
            fprintf(fid13, '%d %d\n', CCSDAT.RCCN(i, j), CCSDAT.ZCCN(i, j));
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
                RGI = CCSDAT.RCCS(i, 2 * I - 1) * F1 + CCSDAT.RCCS(i, 2 * I) * F2 + CCSDAT.RCCS(i, 2 * I + 1) * F3;
                ZGI = CCSDAT.ZCCS(i, 2 * I - 1) * F1 + CCSDAT.ZCCS(i, 2 * I) * F2 + CCSDAT.ZCCS(i, 2 * I + 1) * F3;
                RGI2(count) = RGI;
                ZGI2(count) = ZGI;
                count = count + 1;

                fprintf(fid14, '%d %d\n', RGI, ZGI);
                fprintf(fid15, '%d %d\n', RGI, ZGI);
            end

        end

    end

    CCSDAT.RGI = RGI2;
    CCSDAT.ZGI = ZGI2;

    fclose(fid12);
    fclose(fid13);
    fclose(fid14);
    fclose(fid15);

end

%% maksCCSdata kokomade!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% *************************************************************
%% *************************************************************

%% D_CCS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% D_CCS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% D_CCS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% D_CCS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function CCSDAT = D_CCS(PARAM)

    SR(1:PARAM.CCS, 1) = PARAM.R0 - PARAM.RR;
    SR(1:PARAM.CCS, 3) = SR(1:PARAM.CCS, 1);
    SR(1:PARAM.CCS, 2) = PARAM.R0 + PARAM.RR;
    SZ(1:PARAM.CCS, 1) = PARAM.Z0 + PARAM.RR .* PARAM.CAPPER;
    SZ(1:PARAM.CCS, 2) = PARAM.Z0;
    SZ(1:PARAM.CCS, 3) = PARAM.Z0 - PARAM.RR .* PARAM.CAPPER;

    fid58 = fopen([PARAM.temporary_file_directory '/@D_CCS_SRpoints.txt'], 'w');
    fid59 = fopen([PARAM.temporary_file_directory '/@D_CCS_CheckWrite.txt'], 'w');
    fid13 = fopen([PARAM.temporary_file_directory '/DiscontinuousNodePoints.txt'], 'w');

    for i = 1:PARAM.CCS
        fprintf(fid58, '%d %d\n', PARAM.R0(i), PARAM.Z0(i));

        for j = 1:3
            fprintf(fid58, '%d %d\n', SR(i, j), SZ(i, j));
        end

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

        fprintf('%d%s%d  %d  %d\n', i, '番目CCS: NE/NCCS/NCCN = ', PARAM.NE(i), CCSDAT.NCCS(i), CCSDAT.NCCN(i));

        for j = 1:CCSDAT.NCCS(i) + 1
            fprintf('%d %d %d\n', j, CCSDAT.RCCS(i, j), CCSDAT.ZCCS(i, j));
            fprintf(fid59, '%d %d\n', CCSDAT.RCCS(i, j), CCSDAT.ZCCS(i, j));

            if (II == KADO)
                fprintf(fid59, '%d %d\n', CCSDAT.RCCS(i, j), CCSDAT.ZCCS(i, j));
            end

        end

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

        for j = 1:CCSDAT.NCCN(i)
            fprintf(fid13, '%d %d\n', CCSDAT.RCCN(i, j), CCSDAT.ZCCN(i, j));
        end

    end

    fclose(fid13);
    fclose(fid58);
    fclose(fid59);
end

%% D_CCS kokomade!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
