% 磁場の計算（旧STARB）
function [BZ BR] = cal_B_starb(sen_r, sen_z, coil_r, coil_z, coil_I)
    Mu = 4 * pi * 1.0e-7;
    BZ = 0; BR = 0;
    for k = 1:length(coil_z)
        tmp1 = (coil_z(k) - sen_z).^2 + (coil_r(k) + sen_r).^2;
        kk = 4 * sen_r .* coil_r(k) ./ tmp1;
        [K, E] = ellipke(kk);

        COFA = (sen_r.^2 - coil_r(k)^2 + (coil_z(k) - sen_z).^2) ./ ...
            ((coil_r(k) - sen_r).^2 + (coil_z(k) - sen_z).^2);
        COFB = (coil_r(k)^2 + sen_r.^2 + (coil_z(k) - sen_z).^2) ./ ...
            ((coil_r(k) - sen_r).^2 + (coil_z(k) - sen_z).^2);
        BZ = BZ + coil_I(k) * Mu * (K - COFA .* E) ./ (sqrt(tmp1) * 2 * pi);
        BR = BR + coil_I(k) * Mu * (coil_z(k) - sen_z) ./ sen_r .* (K - COFB .* E) ./ (sqrt(tmp1) * 2 * pi);
    end
end