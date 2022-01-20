
function [GW, HW, GR, GZ, HR, HZ] = cal_integs(rs, zs, rc1, zc1, rc2, zc2, rc3, zc3)
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
    A = rc3 - 2 * rc2 + rc1;
    B = (rc3 - rc1) / 2;
    C = zc3 - 2 * zc2 + zc1;
    D = (zc3 - zc1) / 2;
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
    % rc,zc:積分点
    rc = rc1 * delta1 + rc2 * delta2 + rc3 * delta3;
    zc = zc1 * delta1 + zc2 * delta2 + zc3 * delta3;
    XJA = sqrt((gauss1 * A + B).^2 + (gauss1 * C + D).^2);
    RNOR = (gauss1 * C + D) ./ XJA;
    ZNOR = (-1) * (gauss1 * A + B) ./ XJA;

    [PHI, PHIRs, PHIZs, GSTAR, HSTAR, D] = ...
        cal_starb(rs, zs, rc, zc, RNOR, ZNOR, 1);

    GW(1, :) = ((gauss2 .* XJA .* zeta1) * GSTAR);
    GW(2, :) = ((gauss2 .* XJA .* zeta2) * GSTAR);
    GW(3, :) = ((gauss2 .* XJA .* zeta3) * GSTAR);
    HW(1, :) = ((gauss2 .* XJA .* zeta1) * HSTAR);
    HW(2, :) = ((gauss2 .* XJA .* zeta2) * HSTAR);
    HW(3, :) = ((gauss2 .* XJA .* zeta3) * HSTAR);

    GR(1, :) = ((gauss2 .* XJA .* zeta1) * D.RsG);
    GR(2, :) = ((gauss2 .* XJA .* zeta2) * D.RsG);
    GR(3, :) = ((gauss2 .* XJA .* zeta3) * D.RsG);
    GZ(1, :) = ((gauss2 .* XJA .* zeta1) * D.ZsG);
    GZ(2, :) = ((gauss2 .* XJA .* zeta2) * D.ZsG);
    GZ(3, :) = ((gauss2 .* XJA .* zeta3) * D.ZsG);
    HR(1, :) = ((gauss2 .* XJA .* zeta1) * D.RsH);
    HR(2, :) = ((gauss2 .* XJA .* zeta2) * D.RsH);
    HR(3, :) = ((gauss2 .* XJA .* zeta3) * D.RsH);
    HZ(1, :) = ((gauss2 .* XJA .* zeta1) * D.ZsH);
    HZ(2, :) = ((gauss2 .* XJA .* zeta2) * D.ZsH);
    HZ(3, :) = ((gauss2 .* XJA .* zeta3) * D.ZsH);
end