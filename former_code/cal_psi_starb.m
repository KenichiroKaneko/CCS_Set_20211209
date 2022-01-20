function psi = cal_psi_starb(sen_r, sen_z, coil_r, coil_z, coil_I)
    % 配列を利用した並列化の余地あり
    % 磁束の計算（旧STARB）
    % sen_r,sen_z ... センサーのR位置とZ位置
    % coil_r,coil_z ... 全コイルの 〃
    Mu = 4 * pi * 1e-7;
    psi = 0;
    for k = 1:length(coil_z)
        tmp1 = (coil_z(k) - sen_z).^2 + (coil_r(k) + sen_r).^2;
        kk = 4 .* sen_r .* coil_r(k) ./ tmp1;
        [K, E] = ellipke(kk);

        psi = psi + Mu .* coil_I(k) .* sqrt(sen_r .* coil_r(k)) .* ((1 - kk ./ 2) .* K - E) ./ (pi .* sqrt(kk));
    end
end