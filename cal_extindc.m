function [GW, GR, GZ] = cal_extindc(rs, zs, rc1, zc1, rc2, zc2, rc3, zc3, NONC)
    % rs, zs : センサー位置の配列 1xセンサー数
    % rc, zc : CCS面とは関係ないが他のコードとの整合性で渦電流ノード点位置の変数名をこれにしている
    % NONC = 1 と固定して考えているので、他の場合は考慮していない

    % 通常のGauss積分用のデータ
    gauss1 = [+0.9894009349, -0.9894009349, +0.9445750230, ...
        -0.9445750230, +0.8656312023, -0.8656312023, ...
            +0.7554044083, -0.7554044083, +0.6178762444, ...
            -0.6178762444, +0.4580167776, -0.4580167776, ...
            +0.2816035507, -0.2816035507, +0.0950125098, ...
            -0.0950125098];
    gauss2 = [+0.0271524594, +0.0271524594, +0.0622535239, ...
            +0.0622535239, +0.0951585116, +0.0951585116, ...
            +0.1246289712, +0.1246289712, +0.1495959888, ...
            +0.1495959888, +0.1691565193, +0.1691565193, ...
            +0.1826034150, +0.1826034150, +0.1894506104, ...
            +0.1894506104];

    A = rc3 - 2 * rc2 + rc1;
    B = (rc3 - rc1) / 2;
    C = zc3 - 2 * zc2 + zc1;
    D = (zc3 - zc1) / 2;

    [DSTMIN, GI0, XCQ, YCQ] = MINDST(rs, zs, rc1, zc1, rc2, zc2, rc3, zc3);
    DSTMAX = 0.30;
    whos
    error('error description')
    % ３節点 (rc1,zc1),(rc2,zc2),(rc3,zc3)のうち最も点(rs,zs)に近いのは？
    % その最短距離が DSTMAX より遠いなら、通常のガウス積分に委ねる。
    % 特異性の心配のない境界要素は通常のガウス積分に任せる
    DDDD = 1.0;
    D1 = sqrt((rc1 - rs).^2 + (zc1 - zs).^2);
    DDDD = D1 .* (D1 < DDDD) + DDDD .* (D1 >= DDDD);
    D2 = sqrt((rc2 - rs).^2 + (zc2 - zs).^2);
    DDDD = D2 .* (D2 < DDDD) + DDDD .* (D2 >= DDDD);
    D3 = sqrt((rc3 - rs).^2 + (zc3 - zs).^2);
    DDDD = D3 .* (D3 < DDDD) + DDDD .* (D3 >= DDDD);

    if (DDDD > DSTMAX)
        [GW, GR, GZ] = EXTINDC01(rs, zs, rc1, zc1, rc2, zc2, rc3, zc3);
        return
    end

    F1Q = 0.75 .* GI0 .* (1.5 .* GI0 - 1) ;
    F2Q = (1 - 1.5 .* GI0) .* (1 + 1.5 .* GI0);
    F3Q = 0.75 .* GI0 .* (1.5 .* GI0 + 1);
    % psi_star
    XJAQ = sqrt((GI0 .* A + B).^2 + (GI0 .* C + D).^2);
    R0G0 = (DSTMIN ./ XJAQ).^2;
    UMG = 1 - GI0;
    UPG = 1 + GI0;
    EBM = UMG.^2 + R0G0;
    EBP = UPG.^2 + R0G0;
    EM = UMG .* log(XJAQ .* sqrt(EBM));
    EP = UPG .* log(XJAQ .* sqrt(EBP));
    TANM = atan(UMG .* XJAQ ./ DSTMIN);
    TANP = atan(UPG .* XJAQ ./ DSTMIN);
    EE0 = EM + EP + (TANM + TANP) .* DSTMIN ./ XJAQ - 2;
    EEE = -EE0 .* rs .* XJAQ ./ (2 * pi);
    GWINT1 = F1Q .* EEE;
    GWINT2 = F2Q .* EEE;
    GWINT3 = F3Q .* EEE;
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % dpsi_star/da and dpsi_star/db
    HADAM = (TANM + TANP) ./ (DSTMIN .* XJAQ); %              ! Hadamard singularity
    EBYE2 = (log(EBM) - log(EBP)) ./ (2 .* XJAQ.^2); %  !  自作のもの  E/E**2
    C2 = A ./ 2;
    C1 = B + A .* GI0;
    C0 = A .* GI0 .* GI0 ./ 2 + B .* GI0 + X2 - rs - A .* R0G0 ./ 2; 
    S2 = C ./ 2;
    S1 = D + C .* GI0;
    S0 = C .* GI0 .* GI0 ./ 2 + D .* GI0 + Y2 - zs - C .* R0G0 ./ 2;
    FACT = rs ./ (2 * pi);
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    EER = EE0 ./ (2 * pi) ./ 2;
    FGZ = FACT .* (2 .* S2 ./ XJAQ.^2 + S1 .* EBYE2 + S0 .* HADAM);
    FGR = FACT .* (2 .* C2 ./ XJAQ.^2 + C1 .* EBYE2 + C0 .* HADAM);
    GZINT1 = F1Q .* XJAQ .* FGZ;
    GZINT2 = F2Q .* XJAQ .* FGZ;
    GZINT3 = F3Q .* XJAQ .* FGZ;
    GRINT1 = F1Q .* XJAQ .* (FGR - EER);
    GRINT2 = F2Q .* XJAQ .* (FGR - EER);
    GRINT3 = F3Q .* XJAQ .* (FGR - EER);

    % -1<GI0<1 の時は境界要素を2分割して積分する
    if(abs(GI0) < 1)
        % Compute the values of the shape functions at the integration points
        F1 = GI.* (GI - 1.0) * 0.5;
        F2= 1.0 - GI.^2;
        F3= GI.* (GI+ 1.0) * 0.5;
        % Compute geometrical properties at the integration points
        XCO = rc1 .* F1 + rc2 .* F2 + rc3 .* F3;
        YCO = zc1 .* F1 + zc2 .* F2 + zc3 .* F3;
        XJA = sqrt((GI .* A + B).^2 + (GI .* C + D).^2);
        GJA = XJAQ .* sqrt((GI - GI0).^2 + (DSTMIN ./ XJAQ).^2);
        ETA1 = (GI .* C + D) ./ XJA;
        ETA2 = (-1) .* (GI .* A + B) ./ XJA;
        FNCTL = (-1) .* rs .* log(GJA)./ (pi * 2);
        FNTLR = (-1) .* log(GJA) ./ (pi * 2);
        FNTCR = FACT .* (C2 ./ XJAQ.^2 + C1 .* (GI - GI0) ./ GJA.^2 + C0 ./ GJA.^2);
        FNTCZ = FACT .* (S2 ./ XJAQ.^2 + S1 .* (GI - GI0) ./ GJA.^2 + S0 ./ GJA.^2);

        [PHI, PHIA, PHIB] ;
        
    else
    end

end


function [DSTMIN, GI0, XCQ, YCQ] = MINDST(rs, zs, rc1, zc1, rc2, zc2, rc3, zc3)
    % DSTMINの最小値をなんかしらの手法で探索する
    % 元コードコピペなのでよく理解できていない
    SPAN = 5000;
    G1M = -SPAN;
    G1P = SPAN;
    %
    IMAX = 20;
    CMAX = IMAX;
    GSTAT = G1M;
    GEND = G1P;
    DSTMIN = 1;
    GMAE = 9.99999999999;
    EPS00 = 1e-12;
    GI0 = 0; % ushiki
    XCQ = 0; % ushiki
    YCQ = 0; % ushiki

    for L = 1:100
        DEL = (GEND - GSTAT) ./ CMAX;
        for i = 1:IMAX+1
            CIM1 = i - 1;
            GII = GSTAT + CIM1 .* DEL;
            F1 = GII .* (GII - 1) * 0.5;
            F2 = 1 - GII.^2;
            F3 = GII .* (GII + 1) * 0.5;
            XCO = rc1 .* F1 + rc2 .* F2 + rc3 .* F3;
            YCO = zc1 .* F1 + zc2 .* F2 + zc3 .* F3;
            EDIST = sqrt((XCO - rs).^2 + (YCO - zs).^2);
            GI0 = GII .* (EDIST < DSTMIN) + GI0 .* (EDIST >= DSTMIN);
            XCQ = XCO .* (EDIST < DSTMIN) + XCQ .* (EDIST >= DSTMIN);
            YCQ = YCO .* (EDIST < DSTMIN) + YCQ .* (EDIST >= DSTMIN);
            DSTMIN = EDIST .* (EDIST < DSTMIN) + DSTMIN .* (EDIST >= DSTMIN);
        end
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
    
end