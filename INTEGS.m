%  ***********************************************************************
%++++ SUBROUTINE (No.3)  INTEGS ++++++++++++++++++++++++++++++++++++++++++
%     通常のGauss積分公式を用いて積分計算を行う。
%  ***********************************************************************
function [HW, GW, GR, GZ, HR, HZ] = INTEGS(XS, YS, X1, Y1, X2, Y2, X3, Y3)
    %
    %//// 通常のGauss積分用のデータ ////////////////////////////////////////
    GI = [+0.9894009349D0, -0.9894009349D0, +0.9445750230D0, ...
        -0.9445750230D0, +0.8656312023D0, -0.8656312023D0, ...
            +0.7554044083D0, -0.7554044083D0, +0.6178762444D0, ...
            -0.6178762444D0, +0.4580167776D0, -0.4580167776D0, ...
            +0.2816035507D0, -0.2816035507D0, +0.0950125098D0, ...
            -0.0950125098D0];
    %
    OME = [+0.0271524594D0, +0.0271524594D0, +0.0622535239D0, ...
        +0.0622535239D0, +0.0951585116D0, +0.0951585116D0, ...
            +0.1246289712D0, +0.1246289712D0, +0.1495959888D0, ...
            +0.1495959888D0, +0.1691565193D0, +0.1691565193D0, ...
            +0.1826034150D0, +0.1826034150D0, +0.1894506104D0, ...
            +0.1894506104D0];
    %///////////////////////////////////////////////////////////////////////
    HW = zeros(1, 3);
    GW = zeros(1, 3);
    GR = zeros(1, 3);
    GZ = zeros(1, 3);
    HR = zeros(1, 3);
    HZ = zeros(1, 3);
    %-----------------------------------------------------------------------
    %!% JACOBIANの中身の計算
    %! J=[{(X3-2*X2+X1)ξ+(X3-X1)/2}^2+{(Y3-2*Y2+Y1)ξ+(Y3-Y1)}^2]^0.5
    %!        A              B             C              D
    %
    A = X3 - 2.D0 .* X2 + X1;
    B = (X3 - X1) ./ 2.D0;
    C = Y3 - 2.D0 .* Y2 + Y1;
    D = (Y3 - Y1) ./ 2.D0;
    DELTA1 = zeros(16);
    DELTA2 = zeros(16);
    DELTA3 = zeros(16);
    ZETA1 = zeros(16);
    ZETA2 = zeros(16);
    ZETA3 = zeros(16);
    XCO = zeros(16);
    YCO = zeros(16);
    XJA = zeros(16);
    RNOR = zeros(16);
    ZNOR = zeros(16);
    GSTAR = zeros(16);
    HSTAR = zeros(16);
    DAG = zeros(16);
    DBG = zeros(16);
    DAH = zeros(16);
    DBH = zeros(16);
    %-----------------------------------------------------------------------
    %         ! 形状関数δの計算
    %         ! DELTA1:δ1=ξ(ξ-1)/2
    %         ! DELTA2:δ2=1-ξ^2
    %         ! DELTA3:δ3=ξ(ξ+1)/2
    %
    DELTA1(1:16) = 0.5D0 * GI(1:16) .* (GI(1:16) - 1);
    DELTA2(1:16) = 1 - GI(1:16).^2;
    DELTA3(1:16) = 0.5D0 * GI(1:16) .* (GI(1:16) + 1);
    %-----------------------------------------------------------------------
    %         ! 内挿関数ζの計算
    %         ! ZETA1:ζ1=(3/4)ξ((3/2)ξ-1)
    %         ! ZETA2:ζ2=(1-(3/2)ξ)*(1+(3/2)ξ)
    %         ! ZETA3:ζ3=(3/4)ξ((3/2)ξ+1)
    %
    ZETA1(1:16) = 0.75D0 * GI(1:16) .* (1.5D0 * GI(1:16) - 1);
    ZETA2(1:16) = (1 - 1.5D0 * GI(1:16)) .* (1 + 1.5D0 * GI(1:16));
    ZETA3(1:16) = 0.75D0 * GI(1:16) .* (1.5D0 * GI(1:16) + 1);
    %-----------------------------------------------------------------------
    %         ! XCO,YCO:積分点
    %
    XCO(1:16) = X1 .* DELTA1(1:16) + X2 .* DELTA2(1:16) + X3 .* DELTA3(1:16);
    YCO(1:16) = Y1 .* DELTA1(1:16) + Y2 .* DELTA2(1:16) + Y3 .* DELTA3(1:16);
    XJA(1:16) = sqrt((GI(1:16) .* A + B).^2 + (GI(1:16) .* C + D).^2);
    RNOR(1:16) = (GI(1:16) .* C + D) ./ XJA(1:16);
    ZNOR(1:16) =- (GI(1:16) .* A + B) ./ XJA(1:16);
    %-----------------------------------------------------------------------
    %         ! GW,HWマトリックスの作成
    [PHI, PHIR, PHIZ, PHIA, PHIB, PIRA, PIRB, PIZA, PIZB, GSTAR(1:16), HSTAR(1:16), DAG(1:16), DBG(1:16), DAH(1:16), ...
        DBH(1:16)] = STARB(1, XS, YS, XCO(1:16), YCO(1:16), RNOR(1:16), ZNOR(1:16)); % OK
    %
    GW(1) = sum(GSTAR(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA1(1:16));
    GW(2) = sum(GSTAR(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA2(1:16));
    GW(3) = sum(GSTAR(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA3(1:16));
    %
    HW(1) = sum(HSTAR(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA1(1:16));
    HW(2) = sum(HSTAR(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA2(1:16));
    HW(3) = sum(HSTAR(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA3(1:16));
    %-2003/01/17------------------------------------------------------------
    GR(1) = sum(DAG(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA1(1:16));
    GR(2) = sum(DAG(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA2(1:16));
    GR(3) = sum(DAG(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA3(1:16));
    %
    GZ(1) = sum(DBG(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA1(1:16));
    GZ(2) = sum(DBG(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA2(1:16));
    GZ(3) = sum(DBG(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA3(1:16));
    %
    HR(1) = sum(DAH(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA1(1:16));
    HR(2) = sum(DAH(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA2(1:16));
    HR(3) = sum(DAH(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA3(1:16));
    %
    HZ(1) = sum(DBH(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA1(1:16));
    HZ(2) = sum(DBH(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA2(1:16));
    HZ(3) = sum(DBH(1:16) .* OME(1:16) .* XJA(1:16) .* ZETA3(1:16));
    %-----------------------------------------------------------------------
    %     for I = 1:16
    % %-----------------------------------------------------------------------
    % %         ! 形状関数δの計算
    % %         ! DELTA1:δ1=ξ(ξ-1)/2
    % %         ! DELTA2:δ2=1-ξ^2
    % %         ! DELTA3:δ3=ξ(ξ+1)/2
    % %
    % 	     DELTA1 = 0.5D0*GI(I)*(GI(I)-1);
    % 	     DELTA2 = 1.-GI(I)^2;
    % 	     DELTA3 = 0.5D0*GI(I)*(GI(I)+1);
    % %-----------------------------------------------------------------------
    % %         ! 内挿関数ζの計算
    % %         ! ZETA1:ζ1=(3/4)ξ((3/2)ξ-1)
    % %         ! ZETA2:ζ2=(1-(3/2)ξ)*(1+(3/2)ξ)
    % %         ! ZETA3:ζ3=(3/4)ξ((3/2)ξ+1)
    % %
    % 	     ZETA1 = 0.75D0*GI(I)*(1.5D0*GI(I)-1);
    % 	     ZETA2 = (1-1.5D0*GI(I))*(1+1.5D0*GI(I));
    % 	     ZETA3 = 0.75D0*GI(I)*(1.5D0*GI(I)+1);
    % %-----------------------------------------------------------------------
    % %         ! XCO,YCO:積分点
    % %
    %    	     XCO = X1*DELTA1+X2*DELTA2+X3*DELTA3;
    % 	     YCO = Y1*DELTA1+Y2*DELTA2+Y3*DELTA3;
    % 	     XJA = sqrt((GI(I)*A+B)^2+(GI(I)*C+D)^2);
    % 	     RNOR =  (GI(I)*C+D)/XJA;
    %          ZNOR = -(GI(I)*A+B)/XJA;
    % %-----------------------------------------------------------------------
    % %         ! GW,HWマトリックスの作成
    %          [PHI,PHIR,PHIZ,PHIA,PHIB,PIRA,PIRB,PIZA,PIZB,GSTAR,HSTAR,DAG,DBG,DAH,...
    %          DBH] = STARB(1,XS,YS,XCO,YCO,RNOR,ZNOR); % OK
    % %
    %          GW(1) = GW(1)+GSTAR*OME(I)*XJA*ZETA1;
    % 	     GW(2) = GW(2)+GSTAR*OME(I)*XJA*ZETA2;
    % 	     GW(3) = GW(3)+GSTAR*OME(I)*XJA*ZETA3;
    % %
    %          HW(1) = HW(1)+HSTAR*OME(I)*XJA*ZETA1;
    % 	     HW(2) = HW(2)+HSTAR*OME(I)*XJA*ZETA2;
    % 	     HW(3) = HW(3)+HSTAR*OME(I)*XJA*ZETA3;
    % %-2003/01/17------------------------------------------------------------
    %          GR(1) = GR(1)+DAG*OME(I)*XJA*ZETA1;
    % 	     GR(2) = GR(2)+DAG*OME(I)*XJA*ZETA2;
    % 	     GR(3) = GR(3)+DAG*OME(I)*XJA*ZETA3;
    % %
    % 	     GZ(1) = GZ(1)+DBG*OME(I)*XJA*ZETA1;
    % 	     GZ(2) = GZ(2)+DBG*OME(I)*XJA*ZETA2;
    % 	     GZ(3) = GZ(3)+DBG*OME(I)*XJA*ZETA3;
    % %
    % 	     HR(1) = HR(1)+DAH*OME(I)*XJA*ZETA1;
    % 	     HR(2) = HR(2)+DAH*OME(I)*XJA*ZETA2;
    % 	     HR(3) = HR(3)+DAH*OME(I)*XJA*ZETA3;
    % %
    %    	     HZ(1) = HZ(1)+DBH*OME(I)*XJA*ZETA1;
    % 	     HZ(2) = HZ(2)+DBH*OME(I)*XJA*ZETA2;
    % 	     HZ(3) = HZ(3)+DBH*OME(I)*XJA*ZETA3;
    % %-----------------------------------------------------------------------
end
