%%
function [PHI, PHIR, PHIZ, PHIA, PHIB, PIRA, PIRB, PIZA, PIZB, GSTAR, HSTAR, DAG, DBG, DAH, DBH] = STARB(ID, A, B, R, Z, RNOR, ZNOR)
    P0 = 1.00000000000;
    P1 = 0.44325141463;
    P2 = 0.06260601220;
    P3 = 0.04757383546;
    P4 = 0.01736506451;
    Q1 = 0.24998368310;
    Q2 = 0.09200180037;
    Q3 = 0.04069697526;
    Q4 = 0.00526449639;
    A0 = 1.38629436112;
    A1 = 0.09666344259;
    A2 = 0.03590092383;
    A3 = 0.03742563713;
    A4 = 0.01451196212;
    B0 = 0.50000000000;
    B1 = 0.12498593597;
    B2 = 0.06880248576;
    B3 = 0.03328355346;
    B4 = 0.00441787012;
    DAG = 0; %ushiki
    DBG = 0; %ushiki
    DAH = 0; %ushiki
    DBH = 0; %ushiki
    COFA = 0; %ushiki
    COFB = 0; %ushiki
    PHIA = 0; %ushiki
    PHIB = 0; %ushiki
    BMB1 = 0; %ushiki
    BMB2 = 0; %ushiki
    E1 = 0; %ushiki
    PIRA = 0; %ushiki
    PIRB = 0; %ushiki
    E2 = 0; %ushiki
    PIZA = 0; %ushiki
    PIZB = 0; %ushiki

    % ------------------------------------------------
    %!  Calculate 'k' for the elliptic integral
    RPA = R + A;
    ZMB = Z - B;
    RAB = RPA .* RPA + ZMB .* ZMB;
    AR = A .* R;
    WSK = 4.0D0 .* AR ./ RAB;
    SK = sqrt(WSK);
    Y = 1.0D0 - WSK;
    % ------------------------------------------------
    DLOGY = log(Y);
    %!  Elliptic integral  E(m)
    DE = P0 + Y .* (P1 + Y .* (P2 + Y .* (P3 + Y .* P4))) - DLOGY .* Y .* (Q1 + Y .* (Q2 + Y .* (Q3 + Y .* Q4)));
    %!  Elliptic integral  K(m)
    DK = A0 + Y .* (A1 + Y .* (A2 + Y .* (A3 + Y .* A4))) - DLOGY .* (B0 + Y .* (B1 + Y .* (B2 + Y .* (B3 + Y .* B4))));
    % ------------------------------------------------
    %!     Ψ.* (=PHI)
    PHI = sqrt(AR) .* ((1.0D0 - WSK ./ 2.0D0) .* DK - DE) ./ (pi .* SK);
    %=============
    %! Derivatives of fundamental sol. in the R and Z directions
    RMA = R - A;
    WZMB = ZMB .* ZMB;
    RAZB = RMA .* RMA + WZMB;
    WR = R .* R;
    WA = A .* A;
    COFZ = (WR + WA + WZMB) ./ RAZB;
    COFR = (WR - WA + WZMB) ./ RAZB;
    SRAB = sqrt(RAB);
    BMBO = 2.0D0 .* pi .* SRAB;
    %!     dΨ.*./dr (=PHIR)
    PHIR = R .* (DK - COFR .* DE) ./ BMBO;
    %!     dΨ.*./dz (=PHIZ)
    PHIZ = ZMB .* (DK - COFZ .* DE) ./ BMBO;
    %!============= Normal derivative of the fundamental solution
    HSTAR = (RNOR .* PHIR + ZNOR .* PHIZ) ./ R;
    GSTAR = PHI ./ R;
    %=============
    if (ID > 0)
        COFA = (WA - WR + WZMB) ./ RAZB;
        COFB = COFZ;
        %    !     dΨ.*./da (=PHIA)
        PHIA = A .* (DK - COFA .* DE) ./ BMBO;
        %    !     dΨ.*./db (=PHIB)
        PHIB = -ZMB .* (DK - COFB .* DE) ./ BMBO;
        %
        BMB1 = pi .* SRAB .* RAZB;
        BMB2 = 2.0D0 .* BMB1;
        E1 = (RMA .* DK + RPA .* DE) ./ RAB;
        %    !     d(Ψ.*./dr)./da (=PIRA)
        PIRA = -RPA .* PHIR ./ RAB + ((RPA .* RMA + WZMB) .* E1 - 2.0D0 .* R .* DE + 4.0D0 .* A .* WZMB .* DE ./ RAZB) .* R ./ BMB2;
        %    !     d(Ψ.*./dr)./db (=PIRB)
        PIRB = ZMB .* PHIR ./ RAB + (E1 - 2.0D0 .* RMA .* DE ./ RAZB) .* AR .* ZMB ./ BMB1;
        %
        E2 = (DK + DE) ./ RAB;
        %    !     d(Ψ.*./dz)./da (=PIZA)
        PIZA = -RPA .* PHIZ ./ RAB + ((RMA .* RPA + WZMB) .* E2 - 2.0D0 .* DE - 4.0D0 .* A .* RMA .* DE ./ RAZB) .* R .* ZMB ./ BMB2;
        %    !     d(Ψ.*./dz)./db (=PIZB)
        PIZB = (ZMB ./ RAB - 1.0D0 ./ ZMB) .* PHIZ + (E2 - 2.0D0 .* DE ./ RAZB) .* AR .* WZMB ./ BMB1;
        %
        DAG = PHIA ./ R;
        DBG = PHIB ./ R;
        DAH = (RNOR .* PIRA + ZNOR .* PIZA) ./ R;
        DBH = (RNOR .* PIRB + ZNOR .* PIZB) ./ R;
        %! KORE KORE
    end

    AMU = 1.0D0;
    HSTAR = AMU .* HSTAR;
    GSTAR = AMU .* GSTAR;
    %
    DAG = AMU .* DAG;
    DBG = AMU .* DBG;
    DAH = AMU .* DAH;
    DBH = AMU .* DBH;
end

%
