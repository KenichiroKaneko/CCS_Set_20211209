function [GW, GR, GZ] = EXTINDC01(XP, YP, X1, Y1, X2, Y2, X3, Y3)
    % This subroutine computes the HW and GW matrices
    % that relate a node (XP,YP) with a boundary element
    % using Gauss quadrature
    % RA          = Radius
    % RD1,RD2,RDN = Radius derivatives
    % ETA1, ETA2  = Components of the unit normal to the element
    % XCO,YCO     = Integration point along the element
    % XJA         = Jacobian
    GW = zeros(1, 3);
    GR = zeros(1, 3);
    GZ = zeros(1, 3);
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
    %
    A = X3 - 2.0D0 .* X2 + X1;
    B = (X3 - X1) ./ 2.0D0;
    C = Y3 - 2.0D0 .* Y2 + Y1;
    D = (Y3 - Y1) ./ 2.0D0;
    %
    F1 = zeros(1, 16);
    F2 = zeros(1, 16);
    F3 = zeros(1, 16);
    ZETA1 = zeros(1, 16);
    ZETA2 = zeros(1, 16);
    ZETA3 = zeros(1, 16);
    XCO = zeros(1, 16);
    YCO = zeros(1, 16);
    XJA = zeros(1, 16);
    ETA1 = zeros(1, 16);
    ETA2 = zeros(1, 16);
    %! Compute the values of the shape functions at the integration points
    F1(1:16) = GI(1:16) .* (GI(1:16) - 1.0D0) * 0.5D0;
    F2(1:16) = 1.0D0 - GI(1:16).^2;
    F3(1:16) = GI(1:16) .* (GI(1:16) + 1.0D0) * 0.5D0;
    %
    %-----------------------------------------------------------------------
    %             ! 内挿関数ζの計算
    %             ! ZETA1:ζ1=(3/4)ξ((3/2)ξ-1)
    %             ! ZETA2:ζ2=(1-(3/2)ξ)*(1+(3/2)ξ)
    %             ! ZETA3:ζ3=(3/4)ξ((3/2)ξ+1)
    %
    ZETA1(1:16) = 0.75D0 .* GI(1:16) .* (1.5D0 .* GI(1:16) - 1);
    ZETA2(1:16) = (1 - 1.5D0 .* GI(1:16)) .* (1 + 1.5D0 .* GI(1:16));
    ZETA3(1:16) = 0.75D0 .* GI(1:16) .* (1.5D0 .* GI(1:16) + 1);
    %-----------------------------------------------------------------------
    %
    %    ! Compute geometrical properties at the integration points
    XCO(1:16) = X1 .* F1(1:16) + X2 .* F2(1:16) + X3 .* F3(1:16);
    YCO(1:16) = Y1 .* F1(1:16) + Y2 .* F2(1:16) + Y3 .* F3(1:16);
    XJA(1:16) = sqrt((GI(1:16) .* A + B).^2 + (GI(1:16) .* C + D).^2);
    ETA1(1:16) = (GI(1:16) .* C + D) ./ XJA(1:16);
    ETA2(1:16) = -(GI(1:16) .* A + B) ./ XJA(1:16);
    %
    %    ! Compute GW and HW matrices
    [PHI(1:16), PHIR, PHIZ, PHIA(1:16), PHIB(1:16), PIRA, PIRB, PIZA, PIZB, GSTAR, ...
            HSTAR, DAG, DBG, DAH, DBH] = STARB(1, XP, YP, XCO(1:16), YCO(1:16), ETA1(1:16), ETA2(1:16));
    %
    GW(1) = sum(PHI .* OME .* XJA .* ZETA1);
    GW(2) = sum(PHI .* OME .* XJA .* ZETA2);
    GW(3) = sum(PHI .* OME .* XJA .* ZETA3);
    GR(1) = sum(PHIA .* OME .* XJA .* ZETA1);
    GR(2) = sum(PHIA .* OME .* XJA .* ZETA2);
    GR(3) = sum(PHIA .* OME .* XJA .* ZETA3);
    GZ(1) = sum(PHIB .* OME .* XJA .* ZETA1);
    GZ(2) = sum(PHIB .* OME .* XJA .* ZETA2);
    GZ(3) = sum(PHIB .* OME .* XJA .* ZETA3);
