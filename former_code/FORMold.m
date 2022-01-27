% load("vars_FORM")
load("vars_FORM9")
ECI = ExtCOIL.I .* ExtCOIL.N * 1000;
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

% RCCN,ZCCN : CCSノード点の位置
% NCCN : CCSノード点の数
% RCCS,ZCCS : CCSメッシュ点の位置
% NCCS : CCSメッシュ点の数
RCCN = CCSDAT.RCCN;
ZCCN = CCSDAT.ZCCN;
NCCN = CCSDAT.NCCN;
NCCS = CCSDAT.NCCS;
RCCS = CCSDAT.RCCS;
ZCCS = CCSDAT.ZCCS;

% REV,ZEV : 容器の壁の形
REV = WALL.REV;
ZEV = WALL.ZEV;

% 真空容器上の渦電流節点の設定
% KNE=真空容器の分割境界要素数,  KNM=真空容器上のメッシュ点数,  KNN=真空容器上の節点数
KNE = WALL.KNE;
KNN = WALL.KNN;

% 渦電流接点の位置？　中身空っぽです
RES = WALL.RES;
ZES = WALL.ZES;

% 安定化板上の渦電流節点の設定
% KSE=安定化板の分割境界要素数,  KSN=安定化板上の節点数
KSE = WALL.KSE;
KSN = WALL.KSN;

NONC = PARAM.NONC;
%MXCCS=20;
RMYU0 = 4 * pi * 1e-7;
NE = PARAM.NE;
CCS = PARAM.CCS;
ipconst = PARAM.IPCONST;

RNOR = 0;
ZNOR = 0;

fid99 = fopen([PARAM.temporary_file_directory '/MINDIST.txt'], 'w'); %99
fid100 = fopen([PARAM.temporary_file_directory '/SEKIBUNCHECK.PRI'], 'w'); %100

% **********************************************************************
% if (IT > 1)
% else
% **********************************************************************
%一周で元の位置に戻る
% for I = 1:CCS
%     RCCS(I,NCCS(I)+1) = RCCS(I,1);
%     ZCCS(I,NCCS(I)+1) = ZCCS(I,1);
% end
%    !======================================================================
%    !
%    !    FFFFFF
%    !    FF                             | T-PROBE |
%    !    FFFF                VECTOR  FF=| N-PROBE |
%    !    FF                             |FLUX-LOOP|
%    !    FF                             |   CCS   |
%    !
%ï
%  順問題の解FFを作成！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
%    !FLUX-LOOP
mmm = 0;

for L = 1:SENSOR_FLXLP.NUM
    PSIFLX(L) = 0.0D0;
    %        PSIFLX(L) = GETA;

    [PPSIFLX(1:KCMX), PHIR, PHIZ, PPSIA, PPSIB, PIRA, PIRB, PIZA, PIZB, GSTAR, HSTAR, DAG, DBG, DAH, DBH] = ...
    STARB(0, RS(L), ZS(L), RC(1:KCMX), ZC(1:KCMX), RNOR, ZNOR); % OK

    PSIFLX(L) = PSIFLX(L) + sum(PPSIFLX(1:KCMX) .* ECI(1:KCMX) * RMYU0);
    mmm(L) = PSIFLX(L);
    if CONFIG.DevideFlux
        PSIFLX(L) = PSIFLX(L) / SENSOR_FLXLP.R(L)^2 .* 2;
    end

    FF(L + NAPB) = FF(L + NAPB) - PSIFLX(L); %! 下駄処理を含む
    FC(L + NAPB) = PSIFLX(L); %! コイル電流がつくる磁束
end
% figure()
% plot(FF)
% figure()
% plot(mmm)

%    !T-PROBE & N-PROBE
%    fprintf(WAHAHA,'%s\r\n','T-PROBE & N-PROBE   B  caused by external coils');
mmm = 0;
for L = 1:NAPB
    BR(L) = 0.0D0;
    BZ(L) = 0.0D0;

    [PHI, PHIR, PHIZ, PHIA(1:KCMX), PHIB(1:KCMX), PIRA, PIRB, PIZA, PIZB, GSTAR, HSTAR, DAG, DBG, DAH, DBH] = ...
        STARB(1, RS(L + NFLX), ZS(L + NFLX), RC(1:KCMX), ZC(1:KCMX), RNOR, ZNOR); % OK

    BR(L) = sum((-PHIB(1:KCMX) / RS(L + NFLX)) .* ECI(1:KCMX) * RMYU0);
    BZ(L) = sum((PHIA(1:KCMX) / RS(L + NFLX)) .* ECI(1:KCMX) * RMYU0);
    BBB(L) = BR(L) * cos(TET(L + NFLX)) + BZ(L) * sin(TET(L + NFLX));
    %        fprintf(WAHAHA,'%d %d\r\n',L,BBB);
    FF(L) = FF(L) - BBB(L);
    FC(L) = BBB(L); %! コイル電流寄与
end
%    !CCS
%    fprintf(WAHAHA,'%s\r\n','CCS       PSI  caused by external coils');
for III = 1:CCS

    for L = 1:NCCN(III)
        %            PSIC(L) = GETA;
        PSIC(L) = 0;

        [PPSIC(1:KCMX), PHIR, PHIZ, PPSIA, PPSIB, PIRA, PIRB, PIZA, PIZB, GSTAR, HSTAR, DAG, DBG, DAH, DBH] = ...
            STARB(0, RCCN(III, L), ZCCN(III, L), RC(1:KCMX), ZC(1:KCMX), RNOR, ZNOR); % OK

        PSIC(L) = PSIC(L) + sum(PPSIC(1:KCMX) .* ECI(1:KCMX) * RMYU0);
        mmm(L) = PSIC(L);
        FF(L + sum(NCCN(1:III - 1)) + NAPB + NFLX) =- PSIC(L); %! 下駄処理を含む
        FC(L + sum(NCCN(1:III - 1)) + NAPB + NFLX) = PSIC(L); %! コイル電流寄与
        %            fprintf(WAHAHA,'%d %d %d\r\n',L,PSIC(L),FF(L+sum(NCCN(1:III-1))+NAPB+NFLX));
    end

end

%%
%!======================================================================
%!
%!    AA
%!   AAAA
%!  AA  AA                  |GT -HT|       | T-PROBE |
%! AA    AA       MATRIX AA=|GN -HN|  <--- | N-PROBE |
%! AAAAAAAA                 |GF -HF|       |FLUX-LOOP|
%! AA    AA                 |GC -HC|       |   CCS   |
%! AA    AA
%!
%%%
% fid99 = fopen([PARAM.temporary_file_directory '/MINDIST.txt'], 'w'); %99
% frewind(fid99);
% fid100 = fopen([PARAM.temporary_file_directory '/SEKIBUNCHECK.PRI'], 'w'); %100
% frewind(fid100);
% fprintf(fid99, '%s\n', '****************************************************');
% fprintf(fid99, '%s\n', '***    In the Subr. FORM ***************************');
% fprintf(fid99, '%s\n', '****************************************************');
% fprintf(fid100, '%s\n', '***************************************************');
% fprintf(fid100, '%s\n', '***    In the Subr. FORM **************************');
% fprintf(fid100, '%s\n', '***************************************************');
%%%%
%!FLUX-LOOP
size(AA);
oi = zeros(size(AA));
mmm = 0;
L = 1;
for III = 1:CCS

    for I = 1:NFLX

        for K = 1:NE(III)
            [HW, GW, GR, GZ, HR, HZ] = INTEGS(RS(I), ZS(I), RCCS(III, 2 * K - 1), ZCCS(III, 2 * K - 1), RCCS(III, 2 * K), ZCCS(III, 2 * K), RCCS(III, 2 * K + 1), ZCCS(III, 2 * K + 1)); % OK
            % mmm(I) = GW(1);
            if CONFIG.DevideFlux
                GW = GW / RS(I)^2 * 2;
                HW = HW / RS(I)^2 * 2;
            end

            for JJ = 1:3
                KK = 3 * (K - 1) + JJ;
                KK + 3 * sum(NE(1:III - 1));
                oi(I + NAPB, KK + 3 * sum(NE(1:III - 1))) = 1;

                AA(I + NAPB, KK + 3 * sum(NE(1:III - 1))) = (AA(I + NAPB, KK + 3 * sum(NE(1:III - 1))) + GW(JJ));
                AA(I + NAPB, KK + 3 * sum(NE(1:III - 1)) + sum(NCCN)) = AA(I + NAPB, KK + 3 * sum(NE(1:III - 1)) + sum(NCCN)) - HW(JJ);
            end
        end % 119
    end
    %C
    %!T-PROBE & N-PROBE
    for I = 1:NAPB
        COST = cos(TET(I + NFLX));
        SINT = sin(TET(I + NFLX));

        for K = 1:NE(III)
            [HW, GW, GR, GZ, HR, HZ] = INTEGS(RS(I + NFLX), ZS(I + NFLX), RCCS(III, 2 * K - 1), ...
                ZCCS(III, 2 * K - 1), RCCS(III, 2 * K), ZCCS(III, 2 * K), RCCS(III, 2 * K + 1), ZCCS(III, 2 * K + 1)); % OK
            mmm(I) = HZ(1);
            for JJ = 1:3
                KK = 3 * (K - 1) + JJ;
                G = -COST * GZ(JJ) / RS(I + NFLX) + SINT * GR(JJ) / RS(I + NFLX);
                H = -COST * HZ(JJ) / RS(I + NFLX) + SINT * HR(JJ) / RS(I + NFLX);
                AA(I, KK + 3 * sum(NE(1:III - 1))) = AA(I, KK + 3 * sum(NE(1:III - 1))) + G;

                AA(I, KK + 3 * sum(NE(1:III - 1)) + sum(NCCN)) = AA(I, KK + 3 * sum(NE(1:III - 1)) + sum(NCCN)) - H;
            end

        end %29

    end
    %
    %!CCS
    for I = 1:NCCN(III)
        HII = 0.0;

        for K = 1:NE(III)

            if and((3 * K) >= I, I >= (3 * K - 2))
                %//// 特異値積分 ///////////////////////////////////////////////////////
                if (I == (3 * K - 2)) % 境界内1番目の節点
                    NODO = 1;
                else
                    %                    if(I == (3*K-1))
                    NODO = 2 .* (I == (3 * K - 1)) +3 .* (I ~= (3 * K - 1)); % 境界内2、3番目の節点
                    %                    else
                    %	                    NODO = 3;
                    %                    end
                end

                [GW, HW] = INLOGSA(RCCN(III, I), ZCCN(III, I), RCCS(III, 2 * K - 1), ZCCS(III, 2 * K - 1), RCCS(III, 2 * K), ...
                    ZCCS(III, 2 * K), RCCS(III, 2 * K + 1), ZCCS(III, 2 * K + 1), NODO); % OK
            else
                %//// 通常の積分 ///////////////////////////////////////////////////////
                [HW, GW, GR, GZ, HR, HZ] = INTEGS(RCCN(III, I), ZCCN(III, I), RCCS(III, 2 * K - 1), ZCCS(III, 2 * K - 1), ...
                RCCS(III, 2 * K), ZCCS(III, 2 * K), RCCS(III, 2 * K + 1), ZCCS(III, 2 * K + 1)); % OK
            end

            %  2021/05/30
            % GW = GW / pi / RCCN(III, I)^2;
            % HW = HW / pi / RCCN(III, I)^2;

            for JJ = 1:3
                KK = 3 * (K - 1) + JJ;
                AA(I + sum(NCCN(1:III - 1)) + NAPB + NFLX, KK + 3 * sum(NE(1:III - 1))) = AA(I + sum(NCCN(1:III - 1)) + NAPB + NFLX, KK + 3 * sum(NE(1:III - 1))) + GW(JJ);
                AA(I + sum(NCCN(1:III - 1)) + NAPB + NFLX, KK + 3 * sum(NE(1:III - 1)) + sum(NCCN)) = AA(I + sum(NCCN(1:III - 1)) + NAPB + NFLX, KK + 3 * sum(NE(1:III - 1)) + sum(NCCN)) - HW(JJ);
            end

        end

        for J = 1:NCCN(III)

            if (J == I)
                continue
            end

            HII = HII + AA(I + sum(NCCN(1:III - 1)) + NAPB + NFLX, J + sum(NCCN(1:III - 1)) + sum(NCCN));
        end % 101

        %        if (HII < -0.001)
        HII = HII + 1.0 .* (HII < -0.001);
        %        end
        AA(I + sum(NCCN(1:III - 1)) + NAPB + NFLX, I + sum(NCCN(1:III - 1)) + sum(NCCN)) = -HII;
    end % 140

end


% ??????????????????????????????????????????????????????????????????
%  渦電流   渦電流   渦電流   渦電流   渦電流   渦電流   渦電流
%    　真空容器　　　　　真空容器　　　　　真空容器　　　　真空容器
%     KNE=No. of boundary elements along the vauum vessel
%     KNN=No. of nodes along the vauum vessel (KNN=KNE*2)
%     (REV(),ZEV())=Eddy Current Nodes on the vacuum vessel
% ??????????????????????????????????????????????????????????????????
%****
JJJJ = sum(NCCN) + sum(NCCN);

if (KNE <= 0) %! GOTO 991
else
    AMYU0 = RMYU0 * 1.0D06;
    %  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %    非適合(Non Conforming)渦電流要素の試作  (if NONC=1)
    if (NONC == 0) %GOTO 990
        fprintf('%s\n', '真空容器上の渦電流が flux loop センサーに作るΨ');

        for I = 1:NFLX
            A = RS(I);
            B = ZS(I);

            for K = 1:2:KNN - 1
                fid99 = "a";
                [GW, GR, GZ] = EXTINDC(A, B, REV(K), ZEV(K), REV(K + 1), ZEV(K + 1), REV(K + 2), ...
                    ZEV(K + 2), NONC, fid99, fid100); % OK

                for JJ = 1:3

                    if CONFIG.DevideFlux
                        EE = GW(JJ) * AMYU0 / A^2 * 2; %  2021/05/30
                    else
                        EE = GW(JJ) * AMYU0; % ! flux*AMYU0
                    end

                    if and(K == KNN - 1, JJ == 3)
                        AA(I + NAPB, sum(NCCN) * 2 + 1) = AA(I + NAPB, sum(NCCN) * 2 + 1) + EE; %        ! ### 2
                    else
                        AA(I + NAPB, sum(NCCN) * 2 + K - 1 + JJ) = AA(I + NAPB, sum(NCCN) * 2 + K - 1 + JJ) + EE; %  ! ### 2
                    end

                end

            end

        end

        %
        fprintf('%s\n', '真空容器上の渦電流が 磁場センサーに作るＢ');

        for I = 1:NAPB
            COST = cos(TET(I + NFLX));
            SINT = sin(TET(I + NFLX));
            A = RS(I + NFLX);
            B = ZS(I + NFLX);

            for K = 1:2:KNN - 1
                [GW, GR, GZ] = EXTINDC(A, B, REV(K), ZEV(K), REV(K + 1), ZEV(K + 1), ...
                    REV(K + 2), ZEV(K + 2), NONC, fid99, fid100); % OK

                for JJ = 1:3
                    EE = (-COST * GZ(JJ) + SINT * GR(JJ)) * AMYU0 / A; %! AMYU0/A
                    %%!!   (-COST,SINT) ------ Need to reconfirm --- OK!!
                    if and(K == KNN - 1, JJ == 3)
                        AA(I, sum(NCCN) * 2 + 1) = AA(I, sum(NCCN) * 2 + 1) + EE; %   ! ### 1
                    else
                        AA(I, sum(NCCN) * 2 + K - 1 + JJ) = AA(I, sum(NCCN) * 2 + K - 1 + JJ) + EE; % ! ### 1
                    end

                end

            end

        end

        %
        fprintf('%s\n', '真空容器上の渦電流が CCSに作るΨ');

        for III = 1:CCS

            for I = 1:NCCN(III)
                A = RCCN(III, I);
                B = ZCCN(III, I);

                for K = 1:2:KNN - 1
                    [GW, GR, GZ] = EXTINDC(A, B, REV(K), ZEV(K), REV(K + 1), ZEV(K + 1), ...
                        REV(K + 2), ZEV(K + 2), NONC, fid99, fid100); % OK

                    for JJ = 1:3
                        EE = GW(JJ) * AMYU0; % ! flux*AMYU0
                        % EE = GW(JJ) * AMYU0 / pi / A^2; %  2021/05/30

                        if and(K == KNN - 1, JJ == 3)
                            AA(I + sum(NCCN(1:III - 1)) + NAPB + NFLX, sum(NCCN) * 2 + 1) = AA(I + sum(NCCN(1:III - 1)) + NAPB + NFLX, sum(NCCN) * 2 + 1) + EE; %! ### 3
                        else
                            AA(I + sum(NCCN(1:III - 1)) + NAPB + NFLX, sum(NCCN) * 2 + K - 1 + JJ) = AA(I + sum(NCCN(1:III - 1)) + NAPB + NFLX, sum(NCCN) * 2 + K - 1 + JJ) + EE; % ! ### 3
                        end

                    end

                end

            end

        end
        

        JJJJ = sum(NCCN) + sum(NCCN) + KNN;
        %
    else
        fprintf('%s\n', '真空容器上の渦電流が flux loop センサーに作るΨ');

        for I = 1:NFLX
            A = RS(I);
            B = ZS(I);

            for K = 1:KNE
                [GW, GR, GZ] = EXTINDC(A, B, REV(2 * K - 1), ZEV(2 * K - 1), REV(2 * K), ZEV(2 * K), REV(2 * K + 1), ...
                    ZEV(2 * K + 1), NONC, fid99, fid100);

                for JJ = 1:3
                    KK = 3 * (K - 1) + JJ;

                    if CONFIG.DevideFlux
                        EE = GW(JJ) * AMYU0 / A^2 * 2; %  2021/05/30
                    else
                        EE = GW(JJ) * AMYU0; %      ! flux*AMYU0
                    end

                    AA(I + NAPB, sum(NCCN) * 2 + KK) = AA(I + NAPB, sum(NCCN) * 2 + KK) + EE; %  ! ### 2
                end

            end

        end

        %
        fprintf('%s\n', '真空容器上の渦電流が 磁場センサーに作るＢ');

        for I = 1:NAPB
            COST = cos(TET(I + NFLX));
            SINT = sin(TET(I + NFLX));
            A = RS(I + NFLX);
            B = ZS(I + NFLX);

            for K = 1:KNE
                [GW, GR, GZ] = EXTINDC(A, B, REV(2 * K - 1), ZEV(2 * K - 1), REV(2 * K), ZEV(2 * K), REV(2 * K + 1), ...
                    ZEV(2 * K + 1), NONC, fid99, fid100);

                for JJ = 1:3
                    KK = 3 * (K - 1) + JJ;
                    EE = (-COST * GZ(JJ) + SINT * GR(JJ)) * AMYU0 / A; %  ! AMYU0/A
                    %!!   (-COST,SINT) ------ Need to reconfirm --- OK!!
                    AA(I, sum(NCCN) * 2 + KK) = AA(I, sum(NCCN) * 2 + KK) + EE; %  ! ### 1
                end

            end

        end

        %
        fprintf('%s\n', '真空容器上の渦電流が CCSに作るΨ');

        for III = 1:CCS

            for I = 1:NCCN(III)
                A = RCCN(III, I);
                B = ZCCN(III, I);

                for K = 1:KNE
                    [GW, GR, GZ] = EXTINDC(A, B, REV(2 * K - 1), ZEV(2 * K - 1), REV(2 * K), ZEV(2 * K), REV(2 * K + 1), ...
                        ZEV(2 * K + 1), NONC, fid99, fid100);

                    for JJ = 1:3
                        KK = 3 * (K - 1) + JJ;
                        EE = GW(JJ) * AMYU0; %     ! flux*AMYU0
                        % EE = GW(JJ) * AMYU0 / pi / A^2; %  2021/05/30
                        AA(I + sum(NCCN(1:III - 1)) + NAPB + NFLX, sum(NCCN) * 2 + KK) = AA(I + sum(NCCN(1:III - 1)) + NAPB + NFLX, sum(NCCN) * 2 + KK) + EE; %  ! ### 3
                    end

                end

            end

        end

        JJJJ = sum(NCCN) + sum(NCCN) + KNE * 3; %
    end

figure('Name', 'old')
contour(AA)
figure('Name', 'old')
plot(AA(:, 30))

end

% ??????????????????????????????????????????????????????????????????
%  渦電流   渦電流   渦電流   渦電流   渦電流   渦電流   渦電流
%    　安定化板　　　　　安定化板　　　　　安定化板　　　　安定化板
%CAUTION!! The stabilizer is not closed in the poloidal direction.
%     KSE=No. of boundary elements along the stabilizer
%     KSN=No. of nodes along the stabilizer (KSN=KNE*2+1)
%     (RES(),ZES())=Eddy Current Nodes on the stabilizer
% ??????????????????????????????????????????????????????????????????
if (KSE <= 0) % 991 to 992
else
    AMYU0 = RMYU0 * 1.0D06; % ! NAMUAMUdabutsu  #2
    fprintf('%s\n', '安定化板上の渦電流が flux loop センサーに作るΨ');

    for I = 1:NFLX
        A = RS(I);
        B = ZS(I);

        for K = 1:2:KSN - 2
            [GW, GR, GZ] = EXTINDC(A, B, RES(K), ZES(K), RES(K + 1), ZES(K + 1), RES(K + 2), ...
                ZES(K + 2), NONC, fid99, fid100); % OK

            for JJ = 1:3
                EE = GW(JJ) * AMYU0; %  ! flux*AMYU0
                % EE = GW(JJ) * AMYU0 / pi / A; %  2021/05/30
                AA(I + NAPB, sum(NCCN) * 2 + KNN + K - 1 + JJ) = AA(I + NAPB, sum(NCCN) * 2 + KNN + K - 1 + JJ) + EE; % ! ### 2
            end

        end

    end

    %CC
    fprintf('%s\n', '安定化板上の渦電流が 磁場センサーに作るＢ');

    for I = 1:NAPB
        COST = cos(TET(I + NFLX));
        SINT = sin(TET(I + NFLX));
        A = RS(I + NFLX);
        B = ZS(I + NFLX);

        for K = 1:2:KSN - 2
            [GW, GR, GZ] = EXTINDC(A, B, RES(K), ZES(K), RES(K + 1), ZES(K + 1), RES(K + 2), ...
                ZES(K + 2), NONC, fid99, fid100); % OK

            for JJ = 1:3
                EE = (-COST * GZ(JJ) + SINT * GR(JJ)) * AMYU0 / A; % ! AMYU0/A
                AA(I, sum(NCCN) * 2 + KNN + K - 1 + JJ) = AA(I, sum(NCCN) * 2 + KNN + K - 1 + JJ) + EE; %! ### 1
            end

        end

    end

    %CC
    fprintf('%s\n', '安定化板上の渦電流が CCSに作るΨ');

    for III = 1:CCS

        for I = 1:NCCN(III)
            A = RCCN(III, I);
            B = ZCCN(III, I);

            for K = 1:2:KSN - 2
                [GW, GR, GZ] = EXTINDC(A, B, RES(K), ZES(K), RES(K + 1), ZES(K + 1), RES(K + 2), ...
                    ZES(K + 2), NONC, fid99, fid100); % OK

                for JJ = 1:3
                    % EE = GW(JJ) * AMYU0; %  ! flux*AMYU0
                    % EE = GW(JJ) * AMYU0 / pi / A; %  2021/05/30
                    AA(I + sum(NCCN(1:III - 1)) + NAPB + NFLX, sum(NCCN) * 2 + KNN + K - 1 + JJ) = AA(I + sum(NCCN(1:III - 1)) + NAPB + NFLX, sum(NCCN) * 2 + KNN + K - 1 + JJ) + EE; % ! ### 3
                end

            end

        end

    end

    JJJJ = sum(NCCN) + sum(NCCN) + sum(KNN) + sum(KSN);
    %
    % ??????????????????????????????????????????????????????????????????
    %

    %    end
    % TOTAL IP (USHIKI)
    if (ipconst == 1)

        for III = 1:CCS

            for K = 1:NE(III)
                [HIP] = INTEIP(RCCS(III, 2 * K - 1), ZCCS(III, 2 * K - 1), ...
                    RCCS(III, 2 * K), ZCCS(III, 2 * K), RCCS(III, 2 * K + 1), ZCCS(III, 2 * K + 1)); % OK
                % 2021/05/30
                % HIP = HIP / pi / RCCS(III, 2 * K - 1);

                for JJ = 1:3
                    KK = 3 * (K - 1) + JJ;
                    AA(sum(NCCN) + NAPB + NFLX + 1, KK + 3 * sum(NE(1:III - 1))) = AA(sum(NCCN) + NAPB + NFLX + 1, KK + 3 * sum(NE(1:III - 1))) + HIP(JJ);
                    AA(sum(NCCN) + NAPB + NFLX + 1, KK + 3 * sum(NE(1:III - 1)) + sum(NCCN)) = 0;
                end

            end % 119

        end

    end

end
