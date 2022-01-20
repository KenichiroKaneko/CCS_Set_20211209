function [PHI, PHIRs, PHIZs, GSTAR, HSTAR, D] = cal_starb(rs, zs, rc, zc, RNOR, ZNOR, mode)
    % rs,zs : センサーなどのノード点 (旧A,B)
    % rc,zc : CCSや渦電流ノード点やコイルの位置 (旧R,Z)
    % mode (0 or 1) : どこで計算をやめるか決める
    % SにCが作る磁場磁束などを計算する
    % 行列の計算サイズは[c, s]
    % PHI : µI をかけると磁束になる
    % PHIRs,PHIZs : µI をかけると磁束になる

    D.none = 0;
    num_s = length(rs);
    num_c = length(rc);

    Rs = ones(num_c, 1) * rs;
    Zs = ones(num_c, 1) * zs;
    Rc = rc' * ones(1, num_s);
    Zc = zc' * ones(1, num_s);
    RRNOR = RNOR' * ones(1, num_s);
    ZZNOR = ZNOR' * ones(1, num_s);

    tmp1 = (Rc + Rs).^2 + (Zc - Zs).^2;
    tmp2 = (Rc - Rs).^2 + (Zc - Zs).^2;
    tmp3 = (sqrt(tmp1) .* 2 .* pi);
    tmp4 = pi .* sqrt(tmp1) .* tmp2;
    kk = 4 .* Rs .* Rc ./ tmp1;
    [K, E] = ellipke(kk);
    PHI = sqrt(Rs .* Rc) .* (( 1 - kk ./ 2 ) .* K - E) ./ (pi .* sqrt(kk));

    COFRs = (Rs.^2 - Rc.^2 + (Zc - Zs).^2) ./ tmp2;
    COFZs = (Rc.^2 + Rs.^2 + (Zc - Zs).^2) ./ tmp2;
    COFRc = (Rc.^2 - Rs.^2 + (Zc - Zs).^2) ./ tmp2;
    COFZc = COFZs;

    % dΨ.*./drc (=PHIRs)
    PHIRs = Rs .* (K - COFRs .* E) ./ tmp3;
    % dΨ.*./dzc (=PHIZs)
    PHIZs = (-1) .* (Zc - Zs) .* (K - COFZs .* E) ./ tmp3;
    % dΨ.*./drs (=PHIRc)
    PHIRc = Rc .* (K - COFRc .* E) ./ tmp3;
    % dΨ.*./dzs (=PHIZc)
    PHIZc = (Zc - Zs) .* (K - COFZc .* E) ./ tmp3;


    HSTAR = (RRNOR .* PHIRc + ZZNOR .* PHIZc) ./ Rc;
    GSTAR = PHI ./ Rc; 

    if mode

        % d(Ψ.*./dr)./drs (=PIRRs)
        PIRRs = (-1) .* (Rc + Rs) .* PHIRc ./ tmp1 + (((Rc + Rs) .* (Rc - Rs) + (Zc - Zs).^2) .* ...
            (((Rc - Rs) .* K + (Rc + Rs) .* E) ./ tmp1) - 2 .* Rc .* E + ...
            4 .* Rs .* (Zc - Zs).^2 .* E ./ tmp2) .* Rc ./ (2 .* tmp4);
        % d(Ψ.*./dr)./dzs (=PIRZs)
        PIRZs = (Zc - Zs) .* PHIRc ./ tmp1 + (((Rc - Rs) .* K + (Rc + Rs) .* E) ./ tmp1 - ...
            2 .* (Rc - Rs) .* E ./ tmp2) .* (Rc .* Rs) .* (Zc - Zs) ./ tmp4;
        % d(Ψ.*./dz)./drs (=PIZRS)
        PIZRs = (-1) .* (Rc + Rs) .* PHIZc ./ tmp1 + (((Rc - Rs) .* (Rc + Rs) + (Zc - Zs).^2) .* (K + E) ./ tmp1 - ...
            2 .* E - 4 .* Rs .* (Rc - Rs) .* E ./ tmp2) .* Rc .* (Zc - Zs) ./ (2 .* tmp4);
        % d(Ψ.*./dz)./dzs (=PIZZS)
        PIZZs = ((Zc - Zs) ./ tmp1 - 1 ./ (Zc - Zs)) .* PHIZc + ((K + E) ./ tmp1 - 2 .* E ./ tmp2) .* ...
            (Rc .* Rs) .* (Zc - Zs).^2 ./ tmp4;

        D.RsG = PHIRs ./ Rc;
        D.ZsG = PHIZs ./ Rc;
        D.RsH = (RRNOR .* PIRRs + ZZNOR .* PIZRs) ./ Rc;
        D.ZsH = (RRNOR .* PIRZs + ZZNOR .* PIZZs) ./ Rc;

    end
end