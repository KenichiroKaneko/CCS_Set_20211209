function [FC, BR, BZ, PSIFLX, PSIC, AA, FF] = ...
        FORM(PARAM, CONFIG, FF, ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL)

    NMAX = SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM + sum(CCSDAT.NCCN);
    JMAX = sum(CCSDAT.NCCN) + sum(CCSDAT.NCCN) + sum(WALL.KNN) + sum(WALL.KSN);
    AA = zeros(NMAX, JMAX);

    ECI = ExtCOIL.I .* ExtCOIL.N * 1000;
    KCMX = ExtCOIL.NUM;
    rs = [SENSOR_TPRB.R SENSOR_NPRB.R SENSOR_FLXLP.R];
    zs = [SENSOR_TPRB.Z SENSOR_NPRB.Z SENSOR_FLXLP.Z];

    rc = ExtCOIL.R;
    zc = ExtCOIL.Z;
    TET = [SENSOR_TPRB.TET SENSOR_NPRB.TET SENSOR_FLXLP.TET];
    NAPB = SENSOR_TPRB.NUM + SENSOR_NPRB.NUM;
    NFLX = SENSOR_FLXLP.NUM;
    ind_FL = NAPB + 1:NAPB + NFLX;
    ind_BZ = 1:NAPB;

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

    % 渦電流接点の位置 渦電流があるところに置いて使うやつ
    RES = WALL.RES;
    ZES = WALL.ZES;

    % 安定化板上の渦電流節点の設定
    % KSE=安定化板の分割境界要素数,  KSN=安定化板上の節点数
    KSE = WALL.KSE;
    KSN = WALL.KSN;

    NONC = PARAM.NONC;
    RMYU0 = 4 * pi * 1e-7;
    NE = PARAM.NE;
    CCS = PARAM.CCS;
    ipconst = PARAM.IPCONST;

    RNOR = 0;
    ZNOR = 0;

    % AA • X = FF
    % 順問題の解FFを作成
    %
    % FFFFFF
    % FF                    | T-PROBE |
    % FFFFFF     VECTOR  FF=| N-PROBE |
    % FF                    |FLUX-LOOP|
    % FF                    |   CCS   |

    %% 磁束センサーに外部コイルが作る磁束を差し引く
    [PHI, PHIRs, PHIZs, GSTAR, HSTAR, D] = ... 
        cal_starb(rs(ind_FL), zs(ind_FL), rc, zc, 0, 0, 0);
    psi = RMYU0 * ECI * PHI;
    if CONFIG.DevideFlux
        psi = psi ./ SENSOR_FLXLP.R.^2 .* 2;
    end
    FF(ind_FL) = FF(ind_FL) - psi;
    FC(ind_FL) = psi;
    PSIFLX = psi;

    %% 磁場センサーに外部コイルが作る磁場を差し引く
    % [BZ BR] = cal_B_starb(rs(ind_BZ), zs(ind_BZ), rc, zc, ECI);
    [PHI, PHIRs, PHIZs, GSTAR, HSTAR, D] = ... 
        cal_starb(rs(ind_BZ), zs(ind_BZ), rc, zc, 0, 0, 0);
    BR = RMYU0 * ECI * (-1) * PHIZs ./ rs(ind_BZ);
    BZ = RMYU0 * ECI * PHIRs ./ rs(ind_BZ);
    B = BR .* cos(TET(ind_BZ)) + BZ .* sin(TET(ind_BZ));
    FF(ind_BZ) = FF(ind_BZ) - B;
    FC(ind_BZ) = B;
    
    %% CCS点に外部コイルが作る磁場を差し引く
    PSIC = [];
    for i = 1:CCS
        % 1~9, 10~18のリスト
        ind_CCS_list = [[1:NCCN(1)];[NCCN(1)+1:NCCN(1)*2]] + NAPB + NFLX;
        [PHI, PHIRs, PHIZs, GSTAR, HSTAR, D] = ... 
            cal_starb(RCCN(i, :), ZCCN(i, :), rc, zc, 0, 0, 0);
        
        psi = RMYU0 * ECI * PHI;
        ind_CCS = ind_CCS_list(i,:);
        FF(ind_CCS) = FF(ind_CCS) - psi;
        FC(ind_CCS) = psi;
        PSIC = psi;
    end

    % 関係式の行列 既知数の数ｘ未知数の数の大きさの行列
    %     A
    %    AAA              |GT -HT|       | T-PROBE |
    %   AA AA   MATRIX AA=|GN -HN|  <--- | N-PROBE |
    %  AAAAAAA            |GF -HF|       |FLUX-LOOP|
    % AA     AA           |GC -HC|       |   CCS   |

    % CCSと磁場磁束センサー・自分以外のCCS点との関係式
    for i = 1:CCS
        for j = 1:NE(i)
            % CCSメッシュ点と磁束センサーの関係式
            % [GW HW] = cal_psi_integs( ...
            % rs(ind_FL), zs(ind_FL), ...
            %     RCCS(i, 2 * j - 1), ZCCS(i, 2 * j - 1), ...
            %     RCCS(i, 2 * j), ZCCS(i, 2 * j), ...
            %     RCCS(i, 2 * j + 1), ZCCS(i, 2 * j + 1) ...
            % );
            [GW, HW, GR, GZ, HR, HZ] = cal_integs( ...
                rs(ind_FL), zs(ind_FL), ...
                RCCS(i, 2 * j - 1), ZCCS(i, 2 * j - 1), ...
                RCCS(i, 2 * j), ZCCS(i, 2 * j), ...
                RCCS(i, 2 * j + 1), ZCCS(i, 2 * j + 1) ...
            );
            if CONFIG.DevideFlux
                GW = GW ./ rs(ind_FL).^2 .* 2;
                HW = HW ./ rs(ind_FL).^2 .* 2;
            end
            ind_CCS_FL1 = NAPB + 1:NAPB + NFLX;
            ind_CCS_FL2 = (1 + 3 * (j - 1):3 * j) + 3 * sum(NE(1:i - 1));
            AA(ind_CCS_FL1, ind_CCS_FL2) = GW';
            AA(ind_CCS_FL1, ind_CCS_FL2 + sum(NCCN)) = (-1) .* HW';

            % CCSメッシュ点と磁場センサーの関係式
            [GW, HW, GR, GZ, HR, HZ] = cal_integs( ...
            rs(ind_BZ), zs(ind_BZ), ...
                RCCS(i, 2 * j - 1), ZCCS(i, 2 * j - 1), ...
                RCCS(i, 2 * j), ZCCS(i, 2 * j), ...
                RCCS(i, 2 * j + 1), ZCCS(i, 2 * j + 1) ...
            );
            ind_CCS_BZ1 = 1:NAPB;
            ind_CCS_BZ2 = ind_CCS_FL2;
            G = (-1) .* cos(TET(ind_BZ)) .* GZ ./ rs(ind_BZ) + sin(TET(ind_BZ)) .* GR ./ rs(ind_BZ);
            H = (-1) .* cos(TET(ind_BZ)) .* HZ ./ rs(ind_BZ) + sin(TET(ind_BZ)) .* HR ./ rs(ind_BZ);
            AA(ind_CCS_BZ1, ind_CCS_BZ2) = G';
            AA(ind_CCS_BZ1, ind_CCS_BZ2 + sum(NCCN)) = -H';

        end
        
        % CCSどうしの関係式
        % 配列で計算しようとしたが、CCSノード点の位置によって積分の仕方が異なるので難しくなるからやめた
        for n = 1:NCCN(i)
            for j = 1:NE(i)
                if and(3 * (j - 1) < n, n <= 3 * j)
                    switch rem(n, 3)
                        case 0
                            NODO = 3;
                        case 1
                            NODO = 1;
                        case 2
                            NODO = 2;
                    end
                    % 特異値積分
                    [GW, HW] = INLOGSA( ...
                    RCCN(i, n), ZCCN(i, n), ...
                        RCCS(i, 2 * j - 1), ZCCS(i, 2 * j - 1), ...
                        RCCS(i, 2 * j), ZCCS(i, 2 * j), ...
                        RCCS(i, 2 * j + 1), ZCCS(i, 2 * j + 1), NODO ...
                    );
                else
                    % 通常の積分
                    [GW, HW, GR, GZ, HR, HZ] = cal_integs( ...
                    RCCN(i, n), ZCCN(i, n), ...
                        RCCS(i, 2 * j - 1), ZCCS(i, 2 * j - 1), ...
                        RCCS(i, 2 * j), ZCCS(i, 2 * j), ...
                        RCCS(i, 2 * j + 1), ZCCS(i, 2 * j + 1) ...
                    );
                end
                ind_CCS_CCS1 = n + sum(NCCN(1:i - 1)) + NAPB + NFLX;
                ind_CCS_CCS2 = (1 + 3 * (j - 1):3 * j) + 3 * sum(NE(1:i - 1));
                AA(ind_CCS_CCS1, ind_CCS_CCS2) = GW';
                AA(ind_CCS_CCS1, ind_CCS_CCS2 + sum(NCCN)) = -HW';
            end

            % 元のコードからコピペした謎の処理
            HII = 0;
            for k = 1:NCCN(i)
                if k == n
                    continue
                end
                HII = HII + AA(n + sum(NCCN(1:i - 1)) + NAPB + NFLX, k + sum(NCCN(1:i - 1)) + sum(NCCN));
            end
            HII = HII + 1 .* (HII < -0.001);
            AA(n + sum(NCCN(1:i - 1)) + NAPB + NFLX, n + sum(NCCN(1:i - 1)) + sum(NCCN)) = -HII;
        end

    end

    % 渦電流ノード点と各センサーとの関係式
    % KNE=真空容器の分割境界要素数
    % KNM=真空容器上のメッシュ点数
    % KNN=真空容器上の節点数
    if (KNE <= 0)
    else
        AMYU0 = RMYU0 * 1e06;
        % 非適合(Non Conforming)渦電流要素
        % for i = 1:2:KNN - 1
        % 配列を利用して高速化を今後やる
        %     [GW, GR, GZ] = cal_extindc( ...
        %         rs(ind_FL), zs(ind_FL), ...
        %         REV(i), ZEV(i), REV(i + 1), ZEV(i + 1), REV(i + 2), ZEV(i + 2), NONC);
        % end
        for i = 1:NFLX
            for k = 1:KNE
                % k
                % KNN
                [GW, GR, GZ] = EXTINDC(rs(NAPB + i), zs(NAPB + i),...
                 REV(2*k-1), ZEV(2*k-1), REV(2*k), ZEV(2*k),...
                  REV(2*k + 1), ZEV(2*k + 1), NONC);
                for j = 1:3
                    KK = 3 * (k - 1) + j;
                    if CONFIG.DevideFlux
                        EE = GW(j) * AMYU0 / rs(NAPB + i)^2 * 2;
                    else
                        EE = GW(j) * AMYU0;
                    end
                    AA(i + NAPB, sum(NCCN) * 2 + KK) = EE;
                end
            end
        end

        for i = 1:NAPB
            for k = 1:KNE
                [GW, GR, GZ] = EXTINDC(rs(i), zs(i),...
                 REV(2*k-1), ZEV(2*k-1), REV(2*k), ZEV(2*k),...
                  REV(2*k + 1), ZEV(2*k + 1), NONC);
                for j = 1:3
                    KK = 3 * (k - 1) + j;
                    EE = ((-1)*cos(TET(i))*GZ(j) + sin(TET(i))*GR(j)) * AMYU0/rs(i);
                    AA(i, sum(NCCN)*2+KK) = EE;
                end
            end
        end

        for iii = 1:CCS
            for i = 1:NCCN(iii)
                for k = 1:KNE
                    [GW, GR, GZ] = EXTINDC(RCCN(iii,i), ZCCN(iii,i),...
                    REV(2*k-1), ZEV(2*k-1), REV(2*k), ZEV(2*k),...
                    REV(2*k + 1), ZEV(2*k + 1), NONC);
                    for j = 1:3
                        KK = 3 * (k - 1) + j;
                        EE = GW(j) * AMYU0;
                        AA(i + sum(NCCN(1:iii-1)) + NAPB + NFLX, sum(NCCN) * 2 + KK) = EE;
                    end
                end
            end
        end
    end

end
