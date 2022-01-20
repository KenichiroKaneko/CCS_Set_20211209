
% 磁束の積分を計算
function [GW HW] = cal_psi_integs(sen_r, sen_z, ccs_r1, ccs_z1, ccs_r2, ccs_z2, ccs_r3, ccs_z3)
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

    %% JACOBIANの中身の計算
    % J=[{(X3-2*X2+X1)ξ+(X3-X1)/2}^2+{(Y3-2*Y2+Y1)ξ+(Y3-Y1)}^2]^0.5
    %        A              B             C              D
    A = ccs_r3 - 2 * ccs_r2 + ccs_r1;
    B = (ccs_r3 - ccs_r1) / 2;
    C = ccs_z3 - 2 * ccs_z2 + ccs_z1;
    D = (ccs_z3 - ccs_z1) / 2;
    % 形状関数δの計算
    % DELTA1:δ1=ξ(ξ-1)/2
    % DELTA2:δ2=1-ξ^2
    % DELTA3:δ3=ξ(ξ+1)/2
    delta1 = 0.5 * gauss1 .* (gauss1 - 1);
    delta2 = 1 - gauss1.^2;
    delta3 = 0.5 * gauss1 .* (gauss1 + 1);
    % 内挿関数ζの計算
    % ZETA1:ζ1=(3/4)ξ((3/2)ξ-1)
    % ZETA2:ζ2=(1-(3/2)ξ)*(1+(3/2)ξ)
    % ZETA3:ζ3=(3/4)ξ((3/2)ξ+1)
    zeta1 = 0.75 * gauss1 .* (1.5 * gauss1 - 1);
    zeta2 = (1 - 1.5 * gauss1) .* (1 + 1.5 * gauss1);
    zeta3 = 0.75 * gauss1 .* (1.5 * gauss1 + 1);
    % co_r,co_z:積分点
    co_r = ccs_r1 * delta1 + ccs_r2 * delta2 + ccs_r3 * delta3;
    co_z = ccs_z1 * delta1 + ccs_z2 * delta2 + ccs_z3 * delta3;
    XJA = sqrt((gauss1 * A + B).^2 + (gauss1 * C + D).^2);
    RNOR = (gauss1 * C + D) ./ XJA;
    ZNOR =- (gauss1 * A + B) ./ XJA;

    num = length(sen_r);
    co_zz = co_z' * ones(1, num);
    co_rr = co_r' * ones(1, num);
    sen_zz = ones(16, 1) * sen_z;
    sen_rr = ones(16, 1) * sen_r;
    RRNOR = RNOR' * ones(1, num);
    ZZNOR = ZNOR' * ones(1, num);

    % メインの計算部分
    tmp1 = (co_rr + sen_rr).^2 + (co_zz - sen_zz).^2;
    tmp2 = (co_rr - sen_rr).^2 + (co_zz - sen_zz).^2;
    kk = 4 * sen_rr .* co_rr ./ tmp1;
    [K, E] = ellipke(kk);
    COFR = (co_rr.^2 - sen_rr.^2 + (co_zz - sen_zz).^2) ./ tmp2;
    COFZ = (co_rr.^2 + sen_rr.^2 + (co_zz - sen_zz).^2) ./ tmp2;
    PHIR = co_rr .* (K - COFR .* E) ./ (sqrt(tmp1) * 2 * pi);
    PHIZ = (co_zz - sen_zz) .* (K - COFZ .* E) ./ (sqrt(tmp1) * 2 * pi);
    HSTAR = (RRNOR .* PHIR + ZZNOR .* PHIZ) ./ co_rr;
    GSTAR = sqrt(sen_rr ./ co_rr) .* ((1 - kk ./ 2) .* K - E) ./ (pi .* sqrt(kk));

    GW(1, :) = ((gauss2 .* XJA .* zeta1) * GSTAR);
    GW(2, :) = ((gauss2 .* XJA .* zeta2) * GSTAR);
    GW(3, :) = ((gauss2 .* XJA .* zeta3) * GSTAR);
    HW(1, :) = ((gauss2 .* XJA .* zeta1) * HSTAR);
    HW(2, :) = ((gauss2 .* XJA .* zeta2) * HSTAR);
    HW(3, :) = ((gauss2 .* XJA .* zeta3) * HSTAR);
end