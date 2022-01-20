function [psi, CCR, CCZ, DELGE, RCCS, ZCCS, XPSI, XBBZ, XBBR] = INTER(PARAM, IGOAL, FFOUT, ExtCOIL,...
     SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL)
    
    % % 再構成する磁束のサイズ
    MINR = 10;
    MAXR = 90;
    MINZ = -100;
    MAXZ = 100;
    ICRE = 1;
    JCRE = 2;
    NINT = 0;
    % % cm単位の範囲の指定
    % MINR = 10.815;
    % MAXR = 88.80;
    % MINZ = -99.85;
    % MAXZ = 99.85;
    % % メッシュの間隔delr,delz
    % ICRE = 9.747920133111480e-02;
    % JCRE = 9.827755905511811e-02;
    % NINT = 0;
    % Nz = 2033;
    % Nr = 800;
    % zmin = -9.985000000000001e-01;
    % zmax = 9.985000000000001e-01;
    % rmin = 1.081500000000000e-01;
    % rmax = 0.8880;
    % delr = 9.747920133111480e-04;
    % delz = 9.827755905511811e-04;
    for I = MINR:ICRE:MAXR
        NCOUNT = 0;
        CCR = I / 100.0;
        for J = MINZ:JCRE:MAXZ
            NINT = NINT + 1;
            NCOUNT = NCOUNT + 1;
            CR(NINT) = CCR;
            CZ(NINT) = J / 100.0;
        end
    end

    GETA = 0;
    
    
    % save("vars_inter")
    n100 = 250; % Max. No. of given data ( > NAPB+NFLX+MXCCS) %250
    n50 = 150; % Max. No. of unknowns ( > 2*MXCCS)
    Nedp = 150;
    MXCCS = 20; % MAX! / NUMBER OF ELEMENTS ON THE CCS
    MXINT = 10000; % MAX! / NUMBER OF INTERNAL POINTS

    ECI = ExtCOIL.I .* ExtCOIL.N * 1000;
    ECIGRP = 0;
    KCMX = ExtCOIL.NUM;
    RS = [SENSOR_FLXLP.R SENSOR_TPRB.R SENSOR_NPRB.R];
    ZS = [SENSOR_FLXLP.Z SENSOR_TPRB.Z SENSOR_NPRB.Z];
    RC = ExtCOIL.R;
    ZC = ExtCOIL.Z;
    %ITYPE=[SENSOR_FLXLP.ITYPE SENSOR_TPRB.ITYPE SENSOR_NPRB.ITYPE];
    TET = [SENSOR_FLXLP.TET SENSOR_TPRB.TET SENSOR_NPRB.TET];
    %NTPB=SENSOR_TPRB.NUM;
    %NNPB=SENSOR_NPRB.NUM;
    NAPB = SENSOR_TPRB.NUM + SENSOR_NPRB.NUM;
    NFLX = SENSOR_FLXLP.NUM;
    RCCN = CCSDAT.RCCN;
    ZCCN = CCSDAT.ZCCN;
    NCCN = CCSDAT.NCCN;
    NCCS = CCSDAT.NCCS;
    RCCS = CCSDAT.RCCS;
    ZCCS = CCSDAT.ZCCS;
    REV = WALL.REV;
    ZEV = WALL.ZEV;
    KNE = WALL.KNE;
    KNN = WALL.KNN;
    RES = WALL.RES;
    ZES = WALL.ZES;
    KSE = WALL.KSE;
    KSN = WALL.KSN;
    %Nedp=150;
    NONC = PARAM.NONC;
    %MXCCS=20;
    RMYU0 = 4 * pi * 1e-7;
    NE = PARAM.NE;
    CCS = PARAM.CCS;
    %   ipconst=PARAM.IPCONST;
    AMYU0 = RMYU0 * 1.0D06; %! NAMUAMUdabutsu  #1

    %%
    PSI = zeros(1, MXINT);
    PSIA = zeros(1, MXINT);
    PSIB = zeros(1, MXINT);
    RCCSR = zeros(1, MXCCS + 1);
    ZCCSR = zeros(1, MXCCS + 1);
    FI = zeros(1, MXCCS);
    DFI = zeros(1, MXCCS); %  ! BOUNDARY CONDITION FI:Ψ, DFI:dΨ/dn
    XPSI = zeros(1, n100);
    XBBR = zeros(1, 2000);
    XBBZ = zeros(1, 2000);
    DELGE = 0; % ushiki
    % %
    % fid99 = fopen([PARAM.temporary_file_directory '/MINDIST.txt'], 'w'); %99
    % frewind(fid99);
    % fid100 = fopen([PARAM.temporary_file_directory '/SEKIBUNCHECK.PRI'], 'w'); %100
    fid99 = 0; %99
    fid100 = 0; %100
    %
    for I = 1:CCS
        RCCS(I, NCCS(I) + 1) = RCCS(I, 1);
        ZCCS(I, NCCS(I) + 1) = ZCCS(I, 1);
    end

    %
    DFI(1:sum(NCCN)) = FFOUT(1:sum(NCCN));
    FI(1:sum(NCCN)) = FFOUT(1 + sum(NCCN):sum(NCCN) + sum(NCCN));
    % fprintf('%s %d\n', 'GETA in INTER =', GETA);

    if (IGOAL > 0) %GOTO 999
    else
        %%  磁場センサーに作るＢ'')')
        for L = 1:NAPB
            A = RS(NFLX + L);
            B = ZS(NFLX + L);
            % *******************************************************************
            [PSIL, PSIA, PSIB] = QINTER(A, B, GETA, RCCS, ZCCS, FFOUT, FI, DFI, n100, ...
                n50, RC, ZC, ECIGRP, ECI, KCMX, RCCN, ZCCN, NCCN, KNE, KNN, REV, ZEV, KSE, KSN, RES, ...
                ZES, Nedp, AMYU0, NONC, fid99, fid100, RMYU0, NE, CCS); % OK
            % *******************************************************************
            XBBR(L) = -PSIB / A;
            XBBZ(L) = PSIA / B;
        end

        II = 0;

        for I = 1:NAPB
            XBBB = XBBR(I) * cos(TET(I + NFLX)) + XBBZ(I) * sin(TET(I + NFLX));
        end

        % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %  Flux loop センサーに作るΨ'')')
        DELGE = 0.0D0;
        KNT = 0;

        for L = 1:NFLX
            KNT = KNT + 1;
            A = RS(L);
            B = ZS(L);
            % *******************************************************************
            [XPSI(L), PSIA, PSIB] = QINTER(A, B, GETA, RCCS, ZCCS, FFOUT, FI, ...
                DFI, n100, n50, RC, ZC, ECIGRP, ECI, KCMX, RCCN, ZCCN, NCCN, KNE, ...
                KNN, REV, ZEV, KSE, KSN, RES, ZES, Nedp, AMYU0, NONC, fid99, fid100, RMYU0, NE, CCS); % OK
            % ********************************************************************
        end
    end

    % *****************************************************************
    % *****************************************************************
    %   磁束分布マップの作成
    % *****************************************************************
    % *****************************************************************
    %  Sum of contributions by CCS & external coils
    %!INTER INTER INTER
    %  コイル電流が内点に作るΨ
    for L = 1:NINT
        A = CR(L);
        B = CZ(L);
        % *******************************************************************
        [PSI(L), PSIA, PSIB] = QINTER(A, B, GETA, RCCS, ZCCS, FFOUT, FI, DFI, ...
            n100, n50, RC, ZC, ECIGRP, ECI, KCMX, RCCN, ZCCN, NCCN, KNE, KNN, REV, ZEV, KSE, ...
            KSN, RES, ZES, Nedp, AMYU0, NONC, fid99, fid100, RMYU0, NE, CCS); % OK
        % *******************************************************************
    end


    CCR = unique(CR);
    CCZ = unique(CZ);
    CCR(1) = [];
    psi = reshape(PSI(1:numel(CCR) * numel(CCZ)), numel(CCZ), numel(CCR));

    %              L = 1:NINT; b
end

%% INTER kokomade!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% *************************************************************
%% *************************************************************

%% QINTER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% QINTER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% QINTER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% QINTER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [PSIS, PSISA, PSISB] = QINTER(AS, BS, GETA, RCCS, ZCCS, FFOUT, FI, DFI, n100, ...
        n50, RC, ZC, ECIGRP, ECI, KCMX, RCCN, ZCCN, NCCN, KNE, KNN, REV, ZEV, KSE, KSN, RES, ZES, ...
        Nedp, AMYU0, NONC, fid99, fid100, RMYU0, NE, CCS)
    %   PSIS = GETA;
    PPSI = zeros(1, KCMX);
    PPSIA = zeros(1, KCMX);
    PPSIB = zeros(1, KCMX);
    %   PSISA = 0.0D0;
    %   PSISB = 0.0D0;
    RNOR = 0; % ushiki
    ZNOR = 0; % ushiki
    %
    %% コイル電流が作る磁場の加算
    [PPSI(1:KCMX), PHIR, PHIZ, PPSIA(1:KCMX), PPSIB(1:KCMX), PIRA, PIRB, PIZA, PIZB, GSTAR, HSTAR, DAG, DBG, DAH, DBH] ...
        = STARB(1, AS, BS, RC(1:KCMX), ZC(1:KCMX), RNOR, ZNOR); % OK
    PSIS = GETA + sum(PPSI(1:KCMX) .* ECI(1:KCMX) .* RMYU0);
    PSISA = sum(PPSIA(1:KCMX) .* ECI(1:KCMX) .* RMYU0);
    PSISB = sum(PPSIB(1:KCMX) .* ECI(1:KCMX) .* RMYU0);
    % ??????????????????????????????????????????????????????????????????
    %% 渦電流寄与の加算 (1) ↓
    % ??????????????????????????????????????????????????????????????????
    %    真空容器上の渦電流が作る
    if (KNE > 0)
        %  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        %    非適合(Non Conforming)渦電流要素の試作  (if NONC=1)  #1
        if (NONC == 0) % GOTO 990

            for K = 1:KNE
                [GW, GR, GZ] = EXTINDC(AS, BS, REV(2 * K - 1), ZEV(2 * K - 1), REV(2 * K), ZEV(2 * K), ...
                    REV(2 * K + 1), ZEV(2 * K + 1), NONC, fid99, fid100); % OK

                for JJ = 1:3
                    KJ2 = 2 * K - 2 + JJ;

                    if (KJ2 > KNN)
                        KJ2 = JJ - 2;
                    end

                    PSIS = PSIS + GW(JJ) * FFOUT(sum(NCCN) * 2 + KJ2) * AMYU0;
                    PSISA = PSISA + GR(JJ) * FFOUT(sum(NCCN) * 2 + KJ2) * AMYU0;
                    PSISB = PSISB + GZ(JJ) * FFOUT(sum(NCCN) * 2 + KJ2) * AMYU0;
                end

            end

            %                K = 1:KNE;
            %                [GW(1:3,K),GR(1:3,K),GZ(1:3,K)] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1),REV(2*K),ZEV(2*K),...
            %                REV(2*K+1),ZEV(2*K+1),NONC,fid99,fid100); % OK
            %                JJ = 1:3;
            % 	           KJ2(K,1) = 2.*K-2+1;
            % 	           KJ2(K,2) = 2.*K-2+2;
            % 	           KJ2(K,3) = 2.*K-2+3;
            %                KJ2(K,1) = (1-2).*(KJ2(K,1) > KNN) + KJ2(K,1).*(KJ2(K,1) <= KNN);
            %                KJ2(K,2) = (2-2).*(KJ2(K,1) > KNN) + KJ2(K,1).*(KJ2(K,1) <= KNN);
            %                KJ2(K,3) = (3-2).*(KJ2(K,1) > KNN) + KJ2(K,1).*(KJ2(K,1) <= KNN);
            % 	           PSIS = PSIS + sum(sum(GW(JJ,K).*FFOUT(NCCN*2+KJ2(K,JJ))*AMYU0));
            % 	           PSISA = PSISA + sum(sum(GR(JJ,K).*FFOUT(NCCN*2+KJ2(K,JJ))*AMYU0));
            % 	           PSISB = PSISB + sum(sum(GZ(JJ,K).*FFOUT(NCCN*2+KJ2(K,JJ))*AMYU0));

        else

            for K = 1:KNE
                [GW, GR, GZ] = EXTINDC(AS, BS, REV(2 * K - 1), ZEV(2 * K - 1), REV(2 * K), ZEV(2 * K), ...
                    REV(2 * K + 1), ZEV(2 * K + 1), NONC, fid99, fid100); % OK
                %                if (abs(REV(2*K) -0.695) < 0.001 )
                %                    [GW_1,GR_1,GZ_1] =  EXTINDC(AS,BS,REV(2*K-1)-0.002,ZEV(2*K-1),REV(2*K)-0.002,ZEV(2*K),REV(2*K+1)-0.002,...
                %                    ZEV(2*K+1),NONC,fid99,fid100); % OK
                %                    [GW_2,GR_2,GZ_2] =  EXTINDC(AS,BS,REV(2*K-1)+0.002,ZEV(2*K-1),REV(2*K)+0.002,ZEV(2*K),REV(2*K+1)+0.002,...
                %                    ZEV(2*K+1),NONC,fid99,fid100); % OK
                %                    GW = GW + GW_1 + GW_2;
                %                    GR = GR + GR_1 + GR_2;
                %                    GZ = GZ + GZ_1 + GZ_2;
                %                end
                %                if (abs(abs(ZEV(2*K)) -0.2925) < 0.005 )
                %                    [GW_1,GR_1,GZ_1] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1)+0.00625,REV(2*K),ZEV(2*K)+0.00625,REV(2*K+1),...
                %                    ZEV(2*K+1)+0.00625,NONC,fid99,fid100); % OK
                %                    [GW_2,GR_2,GZ_2] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1)-0.00625,REV(2*K),ZEV(2*K)-0.00625,REV(2*K+1),...
                %                    ZEV(2*K+1)-0.00625,NONC,fid99,fid100); % OK
                %                    GW = GW + GW_1 + GW_2;
                %                    GR = GR + GR_1 + GR_2;
                %                    GZ = GZ + GZ_1 + GZ_2;
                %                end
                %                if (abs(REV(2*K) -0.1183365) < 0.001)
                %                    if(ZEV(2*K) > 0)
                %                        fugou = 1;
                %                    else
                %                        fugou = -1;
                %                    end
                %                    [GW_1,GR_1,GZ_1] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1)+0.042*fugou,REV(2*K),ZEV(2*K)+0.042*fugou,REV(2*K+1),...
                %                    ZEV(2*K+1)+0.042*fugou,NONC,fid99,fid100); % OK
                %                    [GW_2,GR_2,GZ_2] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1)+0.084*fugou,REV(2*K),ZEV(2*K)+0.084*fugou,REV(2*K+1),...
                %                    ZEV(2*K+1)+0.084*fugou,NONC,fid99,fid100); % OK
                %                    [GW_3,GR_3,GZ_3] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1)+0.126*fugou,REV(2*K),ZEV(2*K)+0.126*fugou,REV(2*K+1),...
                %                    ZEV(2*K+1)+0.126*fugou,NONC,fid99,fid100); % OK
                %                    [GW_4,GR_4,GZ_4] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1)+0.168*fugou,REV(2*K),ZEV(2*K)+0.168*fugou,REV(2*K+1),...
                %                    ZEV(2*K+1)+0.168*fugou,NONC,fid99,fid100); % OK
                %                    GW = GW + GW_1 + GW_2 + GW_3 + GW_4;
                %                    GR = GR + GR_1 + GR_2 + GR_3 + GR_4;
                %                    GZ = GZ + GZ_1 + GZ_2 + GZ_3 + GZ_4;
                %                end
                %                if (abs(ZEV(2*K)) < 0.001 )
                %                    GW = GW*0;
                %                    GR = GR*0;
                %                    GZ = GZ*0;
                %                end
                for JJ = 1:3
                    KJ2 = 3 * (K - 1) + JJ;
                    PSIS = PSIS + GW(JJ) * FFOUT(sum(NCCN) * 2 + KJ2) * AMYU0;
                    PSISA = PSISA + GR(JJ) * FFOUT(sum(NCCN) * 2 + KJ2) * AMYU0;
                    PSISB = PSISB + GZ(JJ) * FFOUT(sum(NCCN) * 2 + KJ2) * AMYU0;
                end

            end

            %            K = 1:KNE;
            %            [GW(1:3,K),GR(1:3,K),GZ(1:3,K)] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1),REV(2*K),ZEV(2*K),...
            %            REV(2*K+1),ZEV(2*K+1),NONC,fid99,fid100); % OK
            %            JJ = 1:3;
            %            KJ2(K,1) = 3*(K-1)+1;
            %            KJ2(K,2) = 3*(K-1)+2;
            %            KJ2(K,3) = 3*(K-1)+3;
            %            PSIS  = PSIS + sum(sum(GW(JJ,K).*FFOUT(NCCN*2+KJ2(K,JJ))*AMYU0));
            % 	       PSISA = PSISA+ sum(sum(GR(JJ,K).*FFOUT(NCCN*2+KJ2(K,JJ))*AMYU0));
            % 	       PSISB = PSISB+ sum(sum(GZ(JJ,K).*FFOUT(NCCN*2+KJ2(K,JJ))*AMYU0));
            %

        end

    else
    end

    %
    %    安定化板上の渦電流が作る
    if (KSE > 0)

        for K = 1:KSE
            [GW, GR, GZ] = EXTINDC(AS, BS, RES(2 * K - 1), ZES(2 * K - 1), RES(2 * K), ZES(2 * K), ...
                RES(2 * K + 1), ZES(2 * K + 1), NONC, fid99, fid100); % OK

            for JJ = 1:3
                KJ2 = 2 * K - 2 + JJ;
                PSIS = PSIS + GW(JJ) * FFOUT(sum(NCCN) * 2 + KNN + KJ2) * AMYU0;
                PSISA = PSISA + GR(JJ) * FFOUT(sum(NCCN) * 2 + KNN + KJ2) * AMYU0;
                PSISB = PSISB + GZ(JJ) * FFOUT(sum(NCCN) * 2 + KNN + KJ2) * AMYU0;
            end

        end

    else
    end

    %%%c        write(IPR,*) PSIS,PSISA,PSISB
    % ??????????????????????????????????????????????????????????????????
    for III = 1:CCS

        for K = 1:NE(III)
            [HW, GW, GR, GZ, HR, HZ] = INTEGS(AS, BS, RCCS(III, 2 * K - 1), ZCCS(III, 2 * K - 1), RCCS(III, 2 * K), ...
                ZCCS(III, 2 * K), RCCS(III, 2 * K + 1), ZCCS(III, 2 * K + 1)); % OK

            for JJ = 1:3
                KJ2 = 3 * (K - 1) + JJ + 3 * sum(NE(1:III - 1));
                PSIS = PSIS + DFI(KJ2) * GW(JJ) - FI(KJ2) * HW(JJ);
                PSISA = PSISA + DFI(KJ2) * GR(JJ) - FI(KJ2) * HR(JJ);
                PSISB = PSISB + DFI(KJ2) * GZ(JJ) - FI(KJ2) * HZ(JJ);
            end

        end

    end

    %     K = 1:NE
    %     [HW,GW,GR,GZ,HR,HZ] = INTEGS(AS,BS,RCCS(2*K-1),ZCCS(2*K-1),RCCS(2*K),...
    %     ZCCS(2*K),RCCS(2*K+1),ZCCS(2*K+1)); % OK
    %         for JJ = 1:3
    %             KJ2 = 3*(K-1)+JJ;
    %             PSIS = PSIS + DFI(KJ2)*GW(JJ)-FI(KJ2)*HW(JJ);
    %             PSISA= PSISA+ DFI(KJ2)*GR(JJ)-FI(KJ2)*HR(JJ);
    %             PSISB= PSISB+ DFI(KJ2)*GZ(JJ)-FI(KJ2)*HZ(JJ);
    %         end
    %     end
end

%% QINTER kokomade!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
