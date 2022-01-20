function [a] = cal_inter(PARAM, FFout, ExtCOIL,...
    SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL)
    % むりや
a = 0;
    % % 再構成する磁束のサイズ
    MINR = 10;
    MAXR = 90;
    MINZ = -100;
    MAXZ = 100;
    ICRE = 1;
    JCRE = 2;
    NINT = 0;    

    % 再構成する磁束のサイズ
    delz = 0.02; minz = -1; maxz = 1;
    delr = 0.01; minr = 0.1; maxr = 0.9;

    % 磁束をマッピングする点群
    zs = repmat([minz:delz:maxz], 1, (maxr-minr)/delr+1);
    rs = reshape(ones((maxz-minz)/delz+1, 1) * [minr:delr:maxr], 1, []);
    Nmesh = length(rs); % NINT

    n100 = 250; % Max. No. of given data ( > NAPB+NFLX+MXCCS) %250
    n50 = 150; % Max. No. of unknowns ( > 2*MXCCS)
    Nedp = 150;
    MXCCS = 20; % MAX! / NUMBER OF ELEMENTS ON THE CCS
    MXINT = 10000; % MAX! / NUMBER OF INTERNAL POINTS

    ECI = ExtCOIL.I .* ExtCOIL.N * 1000;
    KCMX = ExtCOIL.NUM;
    Rs = [SENSOR_TPRB.R SENSOR_NPRB.R SENSOR_FLXLP.R];
    Zs = [SENSOR_TPRB.Z SENSOR_NPRB.Z SENSOR_FLXLP.Z];

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
    DFI(1:sum(NCCN)) = FFout(1:sum(NCCN));
    FI(1:sum(NCCN)) = FFout(1 + sum(NCCN):sum(NCCN) + sum(NCCN));

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
    AMYU0 = RMYU0 * 1.0e6;
    NE = PARAM.NE;
    CCS = PARAM.CCS;

    % ?
    ECIGRP = 0;

    [PSI(L), PSIA, PSIB] = QINTER(rs, zs, RCCS, ZCCS, FFout, FI, DFI, ...
            n100, n50, rc, zc, ECIGRP, ECI, KCMX, RCCN, ZCCN, NCCN, KNE, KNN, REV, ZEV, KSE, ...
            KSN, RES, ZES, Nedp, AMYU0, NONC, RMYU0, NE, CCS);
    
end

function [psi, psiA, psiB] = QINTER(rs, zs, RCCS, ZCCS, FFout, FI, DFI, ...
    n100, n50, rc, zc, ECIGRP, ECI, KCMX, RCCN, ZCCN, NCCN, KNE, KNN, REV, ZEV, KSE, ...
    KSN, RES, ZES, Nedp, AMYU0, NONC, RMYU0, NE, CCS)


    RNOR = 0;
    ZNOR = 0; 

    % コイル電流が作る磁場の加算
    [PHI, PHIRs, PHIZs, GSTAR, HSTAR, D] = cal_starb(rs, zs, rc, zc, RNOR, ZNOR, 1);
    psi = RMYU0 * ECI * PHI;
    psiA = RMYU0 * ECI * PHIRs;
    psiB = RMYU0 * ECI * PHIZs;
    whos

    for K = 1:KNE
        
    end
end
