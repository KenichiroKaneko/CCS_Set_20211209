%% EDDYP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% EDDYP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% EDDYP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% EDDYP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function DISF = EDDYP(FFOUT, PARAM, CONFIG, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL)
    DISI = zeros(1, 594);
    EDI = zeros(1, 594);

    %     ECI=ExtCOIL.I.*ExtCOIL.N*1000;
    %     KCMX=ExtCOIL.NUM;
    RS = [SENSOR_FLXLP.R SENSOR_TPRB.R SENSOR_NPRB.R];
    ZS = [SENSOR_FLXLP.Z SENSOR_TPRB.Z SENSOR_NPRB.Z];
    %     RC=ExtCOIL.R;
    %     ZC=ExtCOIL.Z;
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
    REVN = WALL.REVN;
    ZEVN = WALL.ZEVN;
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
    ipconst = PARAM.IPCONST;
    %
    % ###############################################################
    %  真空容器に沿った渦電流分布の出力   (*** #4 ***)
    % ###############################################################
    fid40 = fopen([PARAM.output_file_directory, '/EddyCurrentProfile.txt'], 'w');

    
    % Not yet modified
    II = 0;
    DIS = 0;
    TOTR = 0;
    % DISF = 0;

    for K = 1:KNE
        II = II + 1;
        SEG1 = (REV(2 * K - 1) - REVN(3 * K - 2))^2 + (ZEV(2 * K - 1) - ZEVN(3 * K - 2))^2;
        DIS = DIS + sqrt(SEG1);
        TOTR = TOTR + FFOUT(sum(NCCN) * 2 + II) * sqrt(SEG1);
        fprintf(fid40, '%d %d\r\n', DIS, FFOUT(sum(NCCN) * 2 + II));
        DISF.dis(3*K) = DIS; DISF.jeddy(3*K) = FFOUT(sum(NCCN) * 2 + II);
        II = II + 1;
        SEG2 = (REVN(3 * K - 2) - REVN(3 * K - 1))^2 + (ZEVN(3 * K - 2) - ZEVN(3 * K - 1))^2;
        DIS = DIS + sqrt(SEG2);
        TOTR = TOTR + FFOUT(sum(NCCN) * 2 + II) * sqrt(SEG2);
        fprintf(fid40, '%d %d\r\n', DIS, FFOUT(sum(NCCN) * 2 + II));
        DISF.dis(3*K+1) = DIS; DISF.jeddy(3*K+1) = FFOUT(sum(NCCN) * 2 + II);
        II = II + 1;
        SEG3 = (REVN(3 * K - 1) - REVN(3 * K))^2 + (ZEVN(3 * K - 1) - ZEVN(3 * K))^2;
        DIS = DIS + sqrt(SEG3);
        TOTR = TOTR + FFOUT(sum(NCCN) * 2 + II) * sqrt(SEG3);
        fprintf(fid40, '%d %d\r\n', DIS, FFOUT(sum(NCCN) * 2 + II));
        DISF.dis(3*K+2) = DIS; DISF.jeddy(3*K+2) = FFOUT(sum(NCCN) * 2 + II);
        SEG4 = (REVN(3 * K) - REV(2 * K + 1))^2 + (ZEVN(3 * K) - ZEV(2 * K + 1))^2;
        DIS = DIS + sqrt(SEG4);
    end


    % fprintf('%s %d\r\n', 'Total reconstructed eddy current =', TOTR);

    if CONFIG.DataType == 0
        % 数値解で実行するときはReferenceを読み込む
        fid42 = fopen([PARAM.input_file_directory '/jeddy.txt'], 'r');

    else
        % 実験値のときはエラーがでるのでとりあえずReferenceを読み込む
        fid42 = fopen([PARAM.input_file_directory  '/UTST_numel_2033/jeddy.txt'], 'r');
    end

    fid43 = fopen([PARAM.output_file_directory, '/Reference_EddyCurrentProfile.txt'], 'w');
    JMAX = 594;
    TOTI = 0.0;
    DISMAE = 0.0D0;
    temp = textscan(fid42, '%s', 'Delimiter', '\n');

    for J = 1:JMAX
        temp_m = strsplit(temp{1}{J});
        DISI(J) = str2double(temp_m(1));
        EDI (J) = str2double(temp_m(2));
        TOTI = TOTI + (DISI(J) - DISMAE) * EDI(J);
        DISMAE = DISI(J);
    end

    % fprintf('%s %d\r\n', 'Total refrerence eddy current =', TOTI);
    FACT = TOTR / TOTI;
    J = 1:JMAX;
    EDI(J) = EDI(J) * FACT;
    fprintf(fid43, '%d %d\r\n', horzcat(DISI(J)', EDI(J)'));
    %
    fid41 = fopen([PARAM.output_file_directory, '/EddyCurrent_DiscontinuousProfile.txt'], 'w');
    NSEG = 100;
    ASEG = NSEG;
    DEL = 2.0D0 / ASEG;
    DIS = 0.0D0;
    RMAE = REV(1);
    ZMAE = ZEV(1);

    % DISF = [];
    i = 1;

    for K = 1:KNE

        for I = 1:NSEG + 1
            GII = -1.0D0 + DEL * (I - 1);
            %-----------------------------------------------------------------------
            %C 内挿関数ζの計算
            % ZETA1:ζ1=(3/4)ξ((3/2)ξ-1)
            % ZETA2:ζ2=(1-(3/2)ξ)*(1+(3/2)ξ)
            % ZETA3:ζ3=(3/4)ξ((3/2)ξ+1)
            %
            ZETA1 = 0.75D0 * GII * (1.5D0 * GII - 1.0D0);
            ZETA2 = (1.0D0 - 1.5D0 * GII) * (1.0D0 + 1.5D0 * GII);
            ZETA3 = 0.75D0 * GII * (1.5D0 * GII + 1.0D0);
            R = REVN(3 * K - 2) * ZETA1 + REVN(3 * K - 1) * ZETA2 + REVN(3 * K) * ZETA3;
            Z = ZEVN(3 * K - 2) * ZETA1 + ZEVN(3 * K - 1) * ZETA2 + ZEVN(3 * K) * ZETA3;
            EDDY = FFOUT(sum(NCCN) * 2 + 3 * K - 2) * ZETA1 + FFOUT(sum(NCCN) * 2 + 3 * K - 1) * ZETA2 + FFOUT(sum(NCCN) * 2 + 3 * K) * ZETA3;
            DIS = DIS + sqrt((R - RMAE)^2 + (Z - ZMAE)^2);
            fprintf(fid41, '%d %d\r\n', DIS, EDDY);
            % DISF.dis(i) = DIS; DISF.jeddy(i) = EDDY; i = i + 1;
            RMAE = R;
            ZMAE = Z;
        end

        fprintf(fid41, '%s\r\n', '==');
    end

    fclose('all');

end
