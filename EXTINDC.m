function [GW, GR, GZ] = EXTINDC(XP, YP, X1, Y1, X2, Y2, X3, Y3, NONC, fid99, fid100)
    I3 = 3;
    ION = 1; % Temporary output (1/0)=(Yes/No)
    % This subroutine computes the HW and GW matrices
    % that relate a node (XP(=a),YP(=b)) with a boundary element
    % using Gauss quadrature
    % RA          = Radius
    % RD1,RD2,RDN = Radius derivatives
    % ETA1, ETA2  = Components of the unit normal to the element
    % XCO,YCO     = Integration point along the element
    % XJA         = Jacobian
    GW = zeros(1, I3);
    GR = zeros(1, I3);
    GZ = zeros(1, I3);
    %
    GI = [+0.9894009349d0, -0.9894009349d0, +0.9445750230d0, ...
        -0.9445750230d0, +0.8656312023d0, -0.8656312023d0, ...
            +0.7554044083d0, -0.7554044083d0, +0.6178762444d0, ...
            -0.6178762444d0, +0.4580167776d0, -0.4580167776d0, ...
            +0.2816035507d0, -0.2816035507d0, +0.0950125098d0, ...
            -0.0950125098d0];
    OME = [+0.0271524594d0, +0.0271524594d0, +0.0622535239d0, ...
            +0.0622535239d0, +0.0951585116d0, +0.0951585116d0, ...
            +0.1246289712d0, +0.1246289712d0, +0.1495959888d0, ...
            +0.1495959888d0, +0.1691565193d0, +0.1691565193d0, ...
            +0.1826034150d0, +0.1826034150d0, +0.1894506104d0, ...
            +0.1894506104d0];
    % 強制的に特異性の打ち消しを回避
    KYOSEI = 0;

    if (KYOSEI > 0)

        if (NONC == 0)
            [GW, GR, GZ] = EXTINDC00(XP, YP, X1, Y1, X2, Y2, X3, Y3);
        else
            [GW, GR, GZ] = EXTINDC01(XP, YP, X1, Y1, X2, Y2, X3, Y3);
        end

        return
    else
    end

    %
    PAI2 = pi * 2.0D0;
    A = X3 - 2.0D0 .* X2 + X1;
    B = (X3 - X1) ./ 2.0D0;
    C = Y3 - 2.0D0 .* Y2 + Y1;
    D = (Y3 - Y1) ./ 2.0D0;
    %
    %    Start of Min. Distance Search (↓)
    [DSTMIN, GI0, XCQ, YCQ] = MINDST(XP, YP, X1, Y1, X2, Y2, X3, Y3);
    %    End of Min. Distance Search (↑)
    %    ガウス積分に委ねるか否かの判定（その１↓) -----------------------------------
    %    最適な DSTMAX の選択は後でよく考える。
    %                             !
    %CC     DSTMAX=0.10D0    !  for RELAX                       !
    DSTMAX = 0.30D0; %   !   for UTST これ以上にしないとセンサー信号の再現値が過度にばらつく
    %
    %    !  ３節点 (X1,Y1),(X2,Y2),(X3,Y3)のうち最も点(XP,YP)に近いのは？   !
    %    !  その最短距離が DSTMAX より遠いなら、通常のガウス積分に委ねる。      !
    DDDD = 1.0D20;
    D1 = sqrt((X1 - XP).^2 + (Y1 - YP).^2);
    DDDD = D1 .* (D1 < DDDD) + DDDD .* (D1 >= DDDD);
    D2 = sqrt((X2 - XP).^2 + (Y2 - YP).^2);
    DDDD = D2 .* (D2 < DDDD) + DDDD .* (D2 >= DDDD);
    D3 = sqrt((X3 - XP).^2 + (Y3 - YP).^2);
    if (DDDD > DSTMAX)
        % fprintf(fid100, 'DSTMAX/DDDD=%d %d\r\n', DSTMAX, DDDD);
        % fprintf(fid100, '%s\r\n', 'Farther thanDSTMAX, then GO TO EXTINDC00 *********>');

        if (NONC == 0)
            [GW, GR, GZ] = EXTINDC00(XP, YP, X1, Y1, X2, Y2, X3, Y3);
        else
            [GW, GR, GZ] = EXTINDC01(XP, YP, X1, Y1, X2, Y2, X3, Y3);
        end

        return
    else
    end

    %     ! 特異性の心配のない境界要素は通常のガウス積分に任せる (↑)
    %     ! ガウス積分に委ねるか否かの判定 (その１↑) ----------------------------
    %CC *****
    %CC *****
    %
    F1Q = GI0 .* (GI0 - 1.0D0) .* 0.5D0 .* (NONC == 0) + 0.75D0 .* GI0 .* (1.5D0 .* GI0 - 1) .* (NONC ~= 0);
    F2Q = (1.0D0 - GI0.^2) .* (NONC == 0) + (1 - 1.5D0 .* GI0) .* (1 + 1.5D0 .* GI0) .* (NONC ~= 0);
    F3Q = GI0 .* (GI0 + 1.0D0) .* 0.5D0 .* (NONC == 0) + 0.75D0 .* GI0 .* (1.5D0 .* GI0 + 1) .* (NONC ~= 0);
    % psi_star
    XJAQ = sqrt((GI0 .* A + B).^2 + (GI0 .* C + D).^2);
    R0G0 = (DSTMIN ./ XJAQ).^2;
    UMG = 1.0D0 - GI0;
    UPG = 1.0D0 + GI0;
    EBM = UMG.^2 + R0G0;
    EBP = UPG.^2 + R0G0;
    EM = UMG .* log(XJAQ .* sqrt(EBM));
    EP = UPG .* log(XJAQ .* sqrt(EBP));
    TANM = atan(UMG .* XJAQ ./ DSTMIN);
    TANP = atan(UPG .* XJAQ ./ DSTMIN);
    EE0 = EM + EP + (TANM + TANP) .* DSTMIN ./ XJAQ - 2.0D0; % I3
    EEE = -EE0 .* XP .* XJAQ ./ PAI2; % psi_star*XJAQ
    GWINT1 = F1Q .* EEE; %  !
    GWINT2 = F2Q .* EEE; % !
    GWINT3 = F3Q .* EEE; %!
    %
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % dpsi_star/da and dpsi_star/db
    HADAM = (TANM + TANP) ./ (DSTMIN .* XJAQ); %              ! Hadamard singularity
    EBYE2 = (log(EBM) - log(EBP)) ./ (2.0D0 .* XJAQ.^2); %  !  自作のもの  E/E**2
    C2 = A ./ 2.0D0; %
    C1 = B + A .* GI0; %
    C0 = A .* GI0 .* GI0 ./ 2.0D0 + B .* GI0 + X2 - XP - A .* R0G0 ./ 2.0D0; %
    S2 = C ./ 2.0D0;
    S1 = D + C .* GI0;
    S0 = C .* GI0 .* GI0 ./ 2.0D0 + D .* GI0 + Y2 - YP - C .* R0G0 ./ 2.0D0;
    FACT = XP ./ PAI2;
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    EER = EE0 ./ PAI2 ./ 2.0D0;
    FGZ = FACT .* (2.0D0 .* S2 ./ XJAQ.^2 + S1 .* EBYE2 + S0 .* HADAM);
    FGR = FACT .* (2.0D0 .* C2 ./ XJAQ.^2 + C1 .* EBYE2 + C0 .* HADAM);
    GZINT1 = F1Q .* XJAQ .* FGZ;
    GZINT2 = F2Q .* XJAQ .* FGZ;
    GZINT3 = F3Q .* XJAQ .* FGZ;
    GRINT1 = F1Q .* XJAQ .* (FGR - EER);
    GRINT2 = F2Q .* XJAQ .* (FGR - EER);
    GRINT3 = F3Q .* XJAQ .* (FGR - EER);
    %
    %
    %     !   -1<GI0<1 の時は境界要素を2分割して積分する。
    ISPLT = (abs(GI0) <= 0.9999D0);

    % if (ION > 0)
    %     fprintf(fid100, 'abs(GI0)/ISPLT=%d %d\r\n', abs(GI0), ISPLT);
    % end

    % %
    if (ISPLT == 0)
        F1 = zeros(1, 16);
        F2 = zeros(1, 16);
        F3 = zeros(1, 16);
        XCO = zeros(1, 16);
        YCO = zeros(1, 16);
        XJA = zeros(1, 16);
        GJA = zeros(1, 16);
        ETA1 = zeros(1, 16);
        ETA2 = zeros(1, 16);
        FNCTL = zeros(1, 16); %   !  'XP' should be multiplied
        FNTLR = zeros(1, 16);
        FNTCR = zeros(1, 16);
        FNTCZ = zeros(1, 16);
        FF1 = zeros(1, 16);
        FF2 = zeros(1, 16);
        FF3 = zeros(1, 16);
        PHI = zeros(1, 16);
        PHIA = zeros(1, 16);
        PHIB = zeros(1, 16);
        %
        %             ! Compute the values of the shape functions at the integration points
        F1(1:16) = GI(1:16) .* (GI(1:16) - 1.0D0) * 0.5D0;
        F2(1:16) = 1.0D0 - GI(1:16).^2;
        F3(1:16) = GI(1:16) .* (GI(1:16) + 1.0D0) * 0.5D0;
        %        ! Compute geometrical properties at the integration points
        XCO(1:16) = X1 .* F1(1:16) + X2 .* F2(1:16) + X3 .* F3(1:16);
        YCO(1:16) = Y1 .* F1(1:16) + Y2 .* F2(1:16) + Y3 .* F3(1:16);
        XJA(1:16) = sqrt((GI(1:16) .* A + B).^2 + (GI(1:16) .* C + D).^2);
        GJA(1:16) = XJAQ .* sqrt((GI(1:16) - GI0).^2 + (DSTMIN ./ XJAQ).^2);
        ETA1(1:16) = (GI(1:16) .* C + D) ./ XJA(1:16);
        ETA2(1:16) =- (GI(1:16) .* A + B) ./ XJA(1:16);
        FNCTL(1:16) = -XP .* log(GJA(1:16)) ./ PAI2; %   !  'XP' should be multiplied
        FNTLR(1:16) = -log(GJA(1:16)) ./ PAI2 ./ 2.0D0;
        FNTCR(1:16) = FACT .* (C2 ./ XJAQ.^2 + C1 .* (GI(1:16) - GI0) ./ GJA(1:16).^2 + C0 ./ GJA(1:16).^2);
        FNTCZ(1:16) = FACT .* (S2 ./ XJAQ.^2 + S1 .* (GI(1:16) - GI0) ./ GJA(1:16).^2 + S0 ./ GJA(1:16).^2);
        %
        %             ! Compute GW and HW matrices
        [PHI(1:16), PHIR, PHIZ, PHIA(1:16), PHIB(1:16), PIRA, PIRB, PIZA, PIZB, GSTAR, ...
            HSTAR, DAG, DBG, DAH, DBH] = STARB(1, XP, YP, XCO(1:16), YCO(1:16), ETA1(1:16), ETA2(1:16));
        %
        FF1(1:16) = F1(1:16) .* (NONC == 0) + 0.75D0 .* GI(1:16) .* (1.5D0 .* GI(1:16) - 1) .* (NONC ~= 0);
        FF2(1:16) = F2(1:16) .* (NONC == 0) + (1 - 1.5D0 .* GI(1:16)) .* (1 + 1.5D0 .* GI(1:16)) .* (NONC ~= 0);
        FF3(1:16) = F3(1:16) .* (NONC == 0) + 0.75D0 .* GI(1:16) .* (1.5D0 .* GI(1:16) + 1) .* (NONC ~= 0);
        %
        GW(1) = sum(OME(1:16) .* (FF1(1:16) .* XJA(1:16) .* PHI(1:16) - F1Q .* XJAQ .* FNCTL(1:16)));
        GW(2) = sum(OME(1:16) .* (FF2(1:16) .* XJA(1:16) .* PHI(1:16) - F2Q .* XJAQ .* FNCTL(1:16)));
        GW(3) = sum(OME(1:16) .* (FF3(1:16) .* XJA(1:16) .* PHI(1:16) - F3Q .* XJAQ .* FNCTL(1:16)));

        % if (ION > 0)
        %     fprintf(fid100, 'Not divided.  ISPLT/GI0=%d %d\r\n', ISPLT, GI0);
        % end

        %
        GR(1) = sum(OME(1:16) .* (FF1(1:16) .* XJA(1:16) .* PHIA(1:16) - F1Q .* XJAQ .* (FNTCR(1:16) + FNTLR(1:16))));
        GR(2) = sum(OME(1:16) .* (FF2(1:16) .* XJA(1:16) .* PHIA(1:16) - F2Q .* XJAQ .* (FNTCR(1:16) + FNTLR(1:16))));
        GR(3) = sum(OME(1:16) .* (FF3(1:16) .* XJA(1:16) .* PHIA(1:16) - F3Q .* XJAQ .* (FNTCR(1:16) + FNTLR(1:16))));
        %
        GZ(1) = sum(OME(1:16) .* (FF1(1:16) .* XJA(1:16) .* PHIB(1:16) - F1Q .* XJAQ .* FNTCZ(1:16)));
        GZ(2) = sum(OME(1:16) .* (FF2(1:16) .* XJA(1:16) .* PHIB(1:16) - F2Q .* XJAQ .* FNTCZ(1:16)));
        GZ(3) = sum(OME(1:16) .* (FF3(1:16) .* XJA(1:16) .* PHIB(1:16) - F3Q .* XJAQ .* FNTCZ(1:16)));

    else
        %    ! ここからは境界要素を2分割して、それぞれをガウス積分する。
        IL = zeros(1, 16);

        for L = 1:2
            GZL = -1.0D0 .* (L == 1) + GI0 .* (L ~= 1);
            GZR = GI0 .* (L == 1) + 1.0D0 .* (L ~= 1);
            %
            GIM = (GZL + GZR) ./ 2.0D0;
            F1 = GIM .* (GIM - 1.0D0) .* 0.5D0;
            F2 = 1.0D0 - GIM.^2;
            F3 = GIM .* (GIM + 1.0D0) .* 0.5D0;
            XCM = X1 .* F1 + X2 .* F2 + X3 .* F3;
            YCM = Y1 .* F1 + Y2 .* F2 + Y3 .* F3;

            XX1 = X1 .* (L == 1) + XCQ .* (L ~= 1);
            XX2 = XCM .* (L == 1) + XCM .* (L ~= 1);
            XX3 = XCQ .* (L == 1) + X3 .* (L ~= 1);
            YY1 = Y1 .* (L == 1) + YCQ .* (L ~= 1);
            YY2 = YCM .* (L == 1) + YCM .* (L ~= 1);
            YY3 = YCQ .* (L == 1) + Y3 .* (L ~= 1);
            GI0M = 1.0D0 .* (L == 1) -1.0D0 .* (L ~= 1);
            %
            AA = XX3 - 2.0D0 .* XX2 + XX1;
            BB = (XX3 - XX1) ./ 2.0D0;
            CC = YY3 - 2.0D0 .* YY2 + YY1;
            DD = (YY3 - YY1) ./ 2.0D0;
            %
            F1 = zeros(1, 16);
            F2 = zeros(1, 16);
            F3 = zeros(1, 16);
            XCO = zeros(1, 16);
            YCO = zeros(1, 16);
            XJA = zeros(1, 16);
            GJA = zeros(1, 16);
            ETA1 = zeros(1, 16);
            ETA2 = zeros(1, 16);
            FNCTL = zeros(1, 16); %   !  'XP' should be multiplied
            FNTLR = zeros(1, 16);
            FNTCR = zeros(1, 16);
            FNTCZ = zeros(1, 16);
            FF1 = zeros(1, 16);
            FF2 = zeros(1, 16);
            FF3 = zeros(1, 16);
            PHI = zeros(1, 16);
            PHIA = zeros(1, 16);
            PHIB = zeros(1, 16);
            TRNS = zeros(1, 16);
            GII = zeros(1, 16);
            XJAQM = zeros(1, 16);
            IL(1:32) = 1:32;
            %            ! Compute the values of the shape functions at the integration points
            F1(1:16) = GI(1:16) .* (GI(1:16) - 1.0D0) .* 0.5D0;
            F2(1:16) = 1.0D0 - GI(1:16).^2;
            F3(1:16) = GI(1:16) .* (GI(1:16) + 1.0D0) .* 0.5D0;
            %            ! Compute geometrical properties at the integration points
            XCO(1:16) = XX1 .* F1(1:16) + XX2 .* F2(1:16) + XX3 .* F3(1:16);
            YCO(1:16) = YY1 .* F1(1:16) + YY2 .* F2(1:16) + YY3 .* F3(1:16);
            %
            TRNS(1:16) = (1.0D0 + GI(1:16)) ./ 2.0D0;
            GII(1:16) = GZL + (GZR - GZL) .* TRNS(1:16);
            %
            %            ! Redefine F1, F2 and F3
            F1(1:16) = GII(1:16) .* (GII(1:16) - 1.0D0) .* 0.5D0;
            F2(1:16) = 1.0D0 - GII(1:16).^2;
            F3(1:16) = GII(1:16) .* (GII(1:16) + 1.0D0) .* 0.5D0;
            %
            XJA(1:16) = sqrt((GII(1:16) .* A + B).^2 + (GII(1:16) .* C + D).^2);
            GJA(1:16) = XJAQ .* sqrt((GII(1:16) - GI0).^2 + (DSTMIN ./ XJAQ).^2);
            ETA1(1:16) = (GII(1:16) .* C + D) ./ XJA(1:16); % unit vector on the bounadry
            ETA2(1:16) =- (GII(1:16) .* A + B) ./ XJA(1:16); % unit vector on the bounadry
            FNCTL(1:16) = -XP .* log(GJA(1:16)) ./ PAI2; %   !  'XP' should be multiplied
            FNTLR(1:16) = -log(GJA(1:16)) ./ PAI2 ./ 2.0D0;
            FNTCR(1:16) = FACT .* (C2 ./ XJAQ.^2 + C1 .* (GII(1:16) - GI0) ./ GJA(1:16).^2 + C0 ./ GJA(1:16).^2);
            FNTCZ(1:16) = FACT .* (S2 ./ XJAQ.^2 + S1 .* (GII(1:16) - GI0) ./ GJA(1:16).^2 + S0 ./ GJA(1:16).^2);
            %            !     ここまで、これはこれでよい。
            %
            %            ! Compute GW and HW matrices
            [PHI(1:16), PHIR, PHIZ, PHIA(1:16), PHIB(1:16), PIRA, PIRB, PIZA, PIZB, ...
                GSTAR, HSTAR, DAG, DBG, DAH, DBH] = STARB(1, XP, YP, XCO(1:16), YCO(1:16), ETA1(1:16), ETA2(1:16));
            %
            %            !     以下で規格化をしておく。
            XJAQM(1:16) = sqrt((GI0M .* AA + BB).^2 + (GI0M .* CC + DD).^2);
            XJA(1:16) = sqrt((GI((1:16)) .* AA + BB).^2 + (GI((1:16)) .* CC + DD).^2);
            %
            FF1(1:16) = F1(1:16) .* (NONC == 0) + 0.75D0 .* GII(1:16) .* (1.5D0 .* GII(1:16) - 1) .* (NONC ~= 0);
            FF2(1:16) = F2(1:16) .* (NONC == 0) + (1 - 1.5D0 .* GII(1:16)) .* (1 + 1.5D0 .* GII(1:16)) .* (NONC ~= 0);
            FF3(1:16) = F3(1:16) .* (NONC == 0) + 0.75D0 .* GII(1:16) .* (1.5D0 .* GII(1:16) + 1) .* (NONC ~= 0);
            %
            GW(1) = GW(1) + sum(OME(1:16) .* (FF1(1:16) .* XJA(1:16) .* PHI(1:16) - F1Q .* XJAQM(1:16) .* FNCTL(1:16)));
            GW(2) = GW(2) + sum(OME(1:16) .* (FF2(1:16) .* XJA(1:16) .* PHI(1:16) - F2Q .* XJAQM(1:16) .* FNCTL(1:16)));
            GW(3) = GW(3) + sum(OME(1:16) .* (FF3(1:16) .* XJA(1:16) .* PHI(1:16) - F3Q .* XJAQM(1:16) .* FNCTL(1:16)));

            % if (ION > 0)
            %     fprintf(fid100, 'Divided.  ISPLT/GI0=%d %d\r\n', ISPLT, GI0);
            % end

            %
            GR(1) = GR(1) + sum(OME(1:16) .* (FF1(1:16) .* XJA(1:16) .* PHIA(1:16) - F1Q .* XJAQM(1:16) .* (FNTCR(1:16) + FNTLR(1:16))));
            GR(2) = GR(2) + sum(OME(1:16) .* (FF2(1:16) .* XJA(1:16) .* PHIA(1:16) - F2Q .* XJAQM(1:16) .* (FNTCR(1:16) + FNTLR(1:16))));
            GR(3) = GR(3) + sum(OME(1:16) .* (FF3(1:16) .* XJA(1:16) .* PHIA(1:16) - F3Q .* XJAQM(1:16) .* (FNTCR(1:16) + FNTLR(1:16))));
            GZ(1) = GZ(1) + sum(OME(1:16) .* (FF1(1:16) .* XJA(1:16) .* PHIB(1:16) - F1Q .* XJAQM(1:16) .* FNTCZ(1:16)));
            GZ(2) = GZ(2) + sum(OME(1:16) .* (FF2(1:16) .* XJA(1:16) .* PHIB(1:16) - F2Q .* XJAQM(1:16) .* FNTCZ(1:16)));
            GZ(3) = GZ(3) + sum(OME(1:16) .* (FF3(1:16) .* XJA(1:16) .* PHIB(1:16) - F3Q .* XJAQM(1:16) .* FNTCZ(1:16)));
        end

        %
        Tokubetsu = -1.0E10;
        %!c      Tokubetsu=1.0E10
        %c       IF(Tokubetsu.GT.0.1E0.and.GI0.EQ.0.0D0) THEN
        % if (Tokubetsu > 0.1E0 && GI0 > 0.1D0 && DSTMIN < 3.0E-3)
        %     %             IPRN = 3;
        %     %             [] = PRNG_WRZ(IPRN,ISPLT,GI,I16,GW,GR,GZ,I3,GWINT1,GWINT2,GWINT3,GRINT1,GRINT2,GRINT3,GZINT1,GZINT2,GZINT3,fid100)
        %     fprintf('%s %d %d\r\n', 'GI0/DSTMIN =', GI0, DSTMIN);
        %     fprintf('%s %d %d %d\r\n', 'S2/S1/S0 =', S2, S1, S0);
        %     fprintf('%s\r\n', 'Sorry, we quit ordinary computing in this case.');
        %     %             CALL EXIT(0)
        % else
        % end

        %
    end

    GW(1) = GW(1) + GWINT1;
    GW(2) = GW(2) + GWINT2;
    GW(3) = GW(3) + GWINT3;
    GR(1) = GR(1) + GRINT1;
    GR(2) = GR(2) + GRINT2;
    GR(3) = GR(3) + GRINT3;
    GZ(1) = GZ(1) + GZINT1;
    GZ(2) = GZ(2) + GZINT2;
    GZ(3) = GZ(3) + GZINT3;
end

%% MINDST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [DSTMIN, GI0, XCQ, YCQ] = MINDST(XP, YP, X1, Y1, X2, Y2, X3, Y3)
    %C
    %C Start of Min. Distance Search (↓)
    %c      G1M=-1.0D0
    %c      G1P=1.0D0
    %C      SPAN=1.0D0
    SPAN = 5000.0D0;
    G1M = -SPAN;
    G1P = SPAN;
    %
    IMAX = 20;
    CMAX = IMAX;
    GSTAT = G1M;
    GEND = G1P;
    % CIM1 = zeros(1,IMAX+1);
    % GII  = zeros(1,IMAX+1);
    % F1  = zeros(1,IMAX+1);
    % F2  = zeros(1,IMAX+1);
    % F3  = zeros(1,IMAX+1);
    % XCO  = zeros(1,IMAX+1);
    % YCO  = zeros(1,IMAX+1);
    % EDIST  = zeros(1,IMAX+1);
    % GI0 = zeros(1,IMAX+2);
    % XCQ = zeros(1,IMAX+2);
    % YCQ = zeros(1,IMAX+2);
    % DSTMIN = zeros(1,IMAX+2);
    % DSTMIN(1) = 1.0D20;
    DSTMIN = 1.0D20;
    GMAE = 9.99999999999D30;
    EPS00 = 1.0E-12;
    GI0 = 0; % ushiki
    XCQ = 0; % ushiki
    YCQ = 0; % ushiki

    for L = 1:1000
        DEL = (GEND - GSTAT) ./ CMAX;

        for I = 1:IMAX + 1
            CIM1 = I - 1;
            GII = GSTAT + CIM1 .* DEL;
            F1 = GII .* (GII - 1.0D0) * 0.5D0;
            F2 = 1.0D0 - GII.^2;
            F3 = GII .* (GII + 1.0D0) * 0.5D0;
            XCO = X1 .* F1 + X2 .* F2 + X3 .* F3;
            YCO = Y1 .* F1 + Y2 .* F2 + Y3 .* F3;
            EDIST = sqrt((XCO - XP).^2 + (YCO - YP).^2);
            GI0 = GII .* (EDIST < DSTMIN) + GI0 .* (EDIST >= DSTMIN);
            XCQ = XCO .* (EDIST < DSTMIN) + XCQ .* (EDIST >= DSTMIN);
            YCQ = YCO .* (EDIST < DSTMIN) + YCQ .* (EDIST >= DSTMIN);
            DSTMIN = EDIST .* (EDIST < DSTMIN) + DSTMIN .* (EDIST >= DSTMIN);
        end

        %     I = 1:IMAX+1;
        %     CIM1(I) = I - 1;
        %     GII(I) = GSTAT + CIM1(I).*DEL;
        %     F1(I) = GII(I).*(GII(I)-1.0D0)*0.5D0;
        %     F2(I) = 1.0D0-GII(I).^2;
        %     F3(I) = GII(I).*(GII(I)+1.0D0)*0.5D0;
        %     XCO(I) = X1.*F1(I)+X2.*F2(I)+X3.*F3(I);
        %     YCO(I) = Y1.*F1(I)+Y2.*F2(I)+Y3.*F3(I);
        %     EDIST(I) = sqrt((XCO(I)-XP).^2+(YCO(I)-YP).^2);
        %     for I = 1:IMAX+1
        %         DSTMIN(I+1) = EDIST(I).*(EDIST(I) < DSTMIN(I)) + DSTMIN(I).*(EDIST(I) >= DSTMIN(I));
        %         GI0(I+1) = GII(I).*(EDIST(I) < DSTMIN(I)) + GI0(I).*(EDIST(I) >= DSTMIN(I));
        %         XCQ(I+1) = XCO(I).*(EDIST(I) < DSTMIN(I)) + XCQ(I).*(EDIST(I) >= DSTMIN(I));
        %         YCQ(I+1) = YCO(I).*(EDIST(I) < DSTMIN(I)) + YCQ(I).*(EDIST(I) >= DSTMIN(I));
        %     end
        %     DSTMIN(1) = DSTMIN(end);
        %     GI0(1) = GI0(end);
        %     XCQ(1) = XCQ(end);
        %     YCQ(1) = YCQ(end);
        %
        %     EPS = abs((GI0(end)-GMAE)./GMAE);
        %     if (EPS < EPS00)
        %         break
        %     end
        %     GMAE = GI0(end);
        %     GSTAT = GI0(end)-DEL;
        %     GEND = GI0(end)+DEL;
        %     GSTAT = GI0(end).*(GI0(end) <= G1M) + GSTAT.*(GI0(end) > G1M);
        %     GEND = GI0(end).*(GI0(end) >= G1P) + GEND.*(GI0(end) < G1P);
        % end
        EPS = abs((GI0 - GMAE) ./ GMAE);

        if (EPS < EPS00)
            break
        end

        GMAE = GI0;
        GSTAT = GI0 - DEL;
        GEND = GI0 + DEL;
        GSTAT = GI0 .* (GI0 <= G1M) + GSTAT .* (GI0 > G1M);
        GEND = GI0 .* (GI0 >= G1P) + GEND .* (GI0 < G1P);
    end


    % if (EPS > EPS00)
    %     fprintf('%s%d\r\n', '???? Notice here that EPS=', EPS);
    % end

    %! End of Min. Distance Search (↑)
    % EPS0 = 0.10;

    % if (DSTMIN < EPS0)
    %     fprintf(fid99, 'XP/YP=%d %d\r\n', XP, YP);
    %     fprintf(fid99, '     DSTMIN/GI0=%d %d\r\n', DSTMIN, GI0);
    %     fprintf(fid99, '        XCQ/YCQ=%d %d\r\n\r\n', XCQ, YCQ);
    % end

    % if (DSTMIN(end) < EPS0)
    %     fprintf(fid99,'XP/YP=%d %d\r\n',XP,YP);
    %     fprintf(fid99,'     DSTMIN/GI0=%d %d\r\n',DSTMIN(end),GI0(end));
    %     fprintf(fid99,'        XCQ/YCQ=%d %d\r\n\r\n',XCQ(end),YCQ(end));
    % end
    %
    % if (abs(GI0) > SPAN * 0.999D0)
    %     fprintf('%s%d\r\n', '**** Notice here that GI0=', GI0);
    % end

    % if (abs(GI0(end)) > SPAN*0.999D0)
    %     fprintf('%s%d\r\n','**** Notice here that GI0=',GI0(end));
    % end
    % DSTMIN_i =DSTMIN(end);
    % GI0_i =GI0(end);
    % XCQ_i =XCQ(end);
    % YCQ_i =YCQ(end);
end
