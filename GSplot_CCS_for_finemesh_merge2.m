function GSplot_CCS_for_finemesh_merge2(dirname)
    % GSplot_CCS_for_finemesh_merge2('z2033_r602')みたいに実行する
    % CCS03cはUTST装置全体の解
    % CCS05は（途中に仮想壁を置いた）部分的な解
    % 部分的な解を全体に拡張する際に、CCS03cのメッシュデータが必要
    % fluxは平衡解、flux0はpsiの初期状態
    close all;

    saveflag = "modan"; % origin or modan or 0
    dispFigure = 0;
    realFlag = 0;
    % 保存するdir
    id = extractBetween(dirname, 2, length(dirname));
    save_dir = "input2/UTST_numel_" + id;
    if not(exist(save_dir, 'dir'))
        mkdir(save_dir);
    end

    FONT = 'Times';
    fs = 8;

    % 全体の解のデータ
    data_dir3c = "GS_input/z2033_r602/";
    % data_dir3c = "GS_input/z201_r101/";
    % data_dir3c = "GS_input/z201_r81/";
    vars3c = load(data_dir3c + "vars");
    env3c = vars3c.param;
    flux3c = vars3c.psi0;
    psi_v_3c = flux3c';

    % 容器の大きさに拡張したい部分的な解のデータ
    data_dir = "GS_input/" + dirname + "/";
    vars = load(data_dir + "vars");
    env = vars.param;
    flux = vars.psi; % GSコードで求めた磁束
    flux0 = vars.psi0; % EFコイルが作る磁束
    psi0 = flux' - flux0'; % コイルが作る磁束を差し引いた磁束
    psiorig0 = flux';

    % グラフの領域2020/12/21
    z = linspace(env3c.zmin, env3c.zmax, env3c.Nz);
    r = linspace(env3c.rmin, env3c.rmax, env3c.Nr);

    delz = z(2) - z(1);
    delr = r(2) - r(1);

    if env.Nz == env3c.Nz
        % 拡張しなくて良い

        psi = flux';
        figure()
        subplot(1, 6, 1)
        v = linspace(-20, 20, 21);
        contour(r, z, flipud(psi' * 1000), v, 'r')
        psi222 = 0;
        jt = 0;
        psi_nocoil = psi0;
    else
        % 拡張する

        % 片方に数値解をつくった磁場のグラフ2020/12/21
        A = size(psi0);
        psi000 = zeros(env3c.Nr, env3c.Nz - A(1, 2));
        psi00 = [psi0, psi000];
        figure
        ah = subplot(1, 6, 1);
        v = linspace(-20, 20, 21);
        psiorig0 = [psiorig0, psi000];
        contour(r, z, flipud(psiorig0' * 1000), v, 'r')
        ah = subplot(1, 6, 2);
        contour(r, z, flipud(psi00' * 1000), v, 'r')
        % 片方に数値解をつくった磁場のグラフここまで2020/12/21

        psi = psi_v_3c * 0;
        psi(1:env.Nr, 1:env.Nz) = psi0;
        psi(1:env.Nr, end - env.Nz + 1:end) = psi(1:env.Nr, end - env.Nz + 1:end) + psi0(:, end:-1:1);
        % 反対側にコピーした磁場のグラフ2020/12/21
        ah = subplot(1, 6, 3);
        contour(r, z, psi' * 1000, v, 'r')
        % 外部コイルが作る磁束を足した磁束
        psi_nocoil = psi;
        psi = psi + psi_v_3c;
        ah = subplot(1, 6, 4);
        contour(r, z, psi' * 1000, v, 'r')
        %

        % error('error description', A1)
        %% Calc virtual eddy current on center
        flux = vars.mu_jt;
        jt = psi_v_3c * 0;
        jt(1:env.Nr, 1:env.Nz) = jt(1:env.Nr, 1:env.Nz) + flux(1:env.Nz, 1:env.Nr)' / (4 * pi * 1e-7);
        jt(1:env.Nr, end:-1:end - env.Nz + 1) = jt(1:env.Nr, end:-1:end - env.Nz + 1) + ...
            flux(1:env.Nz, 1:env.Nr)' / (4 * pi * 1e-7);
        jt_center = -squeeze(psi0(2:env.Nr - 1, env.Nz - 1)) / delz ./ r(2:env.Nr - 1)' / 2 / pi / 4 / pi / 1e-7 * delr;
        size(jt_center)
        [rr zz] = meshgrid(r, z(1:ceil(length(z) / 2)));
        psi_virtualj = rr * 0;

        const = 2 * pi * 4 * pi * 1e-7;

        % error('error description', A1)

        for k = 1:length(jt_center)
            kk1 = 4 * rr * r(k + 1) ./ ((rr + r(k + 1) + 1e-6).^2 + (zz - z(env.Nz - 1)).^2);
            [K1, E1] = ellipke(kk1);
            kk2 = 4 * rr * r(k + 1) ./ ((rr + r(k + 1) + 1e-6).^2 + (zz + z(env.Nz - 1)).^2);
            [K2, E2] = ellipke(kk2);
            denominator = ((1 ./ (pi * sqrt(kk1)) .* ((1 - kk1 / 2) .* K1 - E1)) + (1 ./ (pi * sqrt(kk2)) .* ((1 - kk2 / 2) .* K2 - E2)));
            psi_virtualj = psi_virtualj + const * jt_center(k) .* sqrt(rr * r(k + 1)) .* denominator;
        end

        psi222 = psi_v_3c * 0;
        psi222(1:env3c.Nr, 1:(env3c.Nz + 1) / 2) = psi_virtualj(1:(env3c.Nz + 1) / 2, 1:env3c.Nr)';
        psi222(1:env3c.Nr, end:-1:end - (env3c.Nz + 1) / 2 + 1 + 1) = psi_virtualj(1:(env3c.Nz + 1) / 2 - 1, 1:env3c.Nr)';

        % 仮想壁に流れる電流を差し引く
        % 平衡解2020/12/21
        psi = psi - psi222;
        ah = subplot(1, 6, 5);
        contour(r, z, psi' * 1000, v, 'r')
    end

    % %% UTST凸部の渦電流を差し引く
    % % UTST凸部のZ座標（zUが上zLが下）
    % [M zU_mesh] = min(abs(z - env3c.z3));
    % [M zL_mesh] = min(abs(z - env3c.z2));
    % % 仮想壁部の１つ内側の電流密度
    % jt_wall = -squeeze(psi_nocoil(env.Nr - 1, zL_mesh + 1:zU_mesh - 1)) / delr / r(env.Nr) / 2 / pi / 4 / pi / 1e-7 * delz;
    % jt_wall = jt_wall;
    % [rr zz] = meshgrid(r, z(1:ceil(length(z) / 2)));
    % psi_virtualj = rr * 0;
    % const = 2 * pi * 4 * pi * 1e-7;
    % i = 1;

    % for k = zL_mesh + 1:zU_mesh - 1
    %     kk1 = 4 * rr * r(env.Nr - 1) ./ ((rr + r(env.Nr - 1) + 1e-6).^2 + (zz - z(k)).^2);
    %     [K1, E1] = ellipke(kk1);
    %     denominator = (1 ./ (pi * sqrt(kk1)) .* ((1 - kk1 / 2) .* K1 - E1));
    %     psi_virtualj = psi_virtualj + const * jt_wall(i) .* sqrt(rr * r(env.Nr - 1)) .* denominator;
    %     i = i + 1;
    % end
    

    psi2222 = psi_v_3c * 0;
    % psi2222(1:env3c.Nr, 1:(env3c.Nz + 1) / 2) = psi_virtualj(1:(env3c.Nz + 1) / 2, 1:env3c.Nr)';
    % psi2222(1:env3c.Nr, end:-1:end - (env3c.Nz + 1) / 2 + 1 + 1) = psi_virtualj(1:(env3c.Nz + 1) / 2 - 1, 1:env3c.Nr)';
    % psi = psi - psi2222;
    % ah = subplot(1, 6, 6);
    % contour(r, z, psi' * 1000, v, 'r')
    % %% UTST凸部の渦電流を差し引く計算ここまで

    %% Calc Bz and Br from psi
    zdiff = z(2:end - 1);
    rdiff = r(2:end - 1);
    Bz = zeros(length(r), length(z));
    Br = zeros(length(r), length(z));

    for j = 1:length(z)
        Bz(2:end - 1, j) = (psi(3:end, j) - psi(1:end - 2, j)) ./ (r(3:end) - r(1:end - 2))' / 2 / pi ./ r(2:end - 1)';
    end

    for j = 1:length(r)
        Br(j, 2:end - 1) =- (psi(j, 3:end) - psi(j, 1:end - 2)) ./ (z(3:end) - z(1:end - 2)) / 2 / pi / r(j);
    end

    %% CCSdata
    datanum = 0;

    % 602*2033
    % 1017 ... Nzの中央
    % 1016 ...  Nzの中央-1
    % 2028 ... Nzから5点内側
    % 6    ... Nzから５点内側

    % 2021/04/28 UTST装置のセンサー配置
    % save("psiGSplotTest", "psi");

    % 2021/07/19
    % フレキシブルな座標
    % I...何点内側か
    I = 5;
    zCenter = ceil(env3c.Nz / 2)
    rSharp = ceil((env3c.rmirror - env3c.rmin) / delr) - I
    zSharp = env3c.Nz - I - floor((env3c.zmax - env3c.z3) / delz)

    count = 0;

    for i = zCenter:env3c.Nz - I
        % １中央から下までインボード側
        count = count + 1;
        r_CCS(count) = r(I + 1);
        z_CCS(count) = z(i);
        psi_CCS(count) = psi(I + 1, i);
        Bz_CCS(count) = Bz(I + 1, i);
        Br_CCS(count) = Br(I + 1, i);
    end

    for i = I + 2:rSharp - 1
        % ２底面インボード側からアウトボード側
        count = count + 1;
        r_CCS(count) = r(i);
        z_CCS(count) = z(env3c.Nz - I);
        psi_CCS(count) = psi(i, env3c.Nz - I);
        Bz_CCS(count) = Bz(i, env3c.Nz - I);
        Br_CCS(count) = Br(i, env3c.Nz - I);
    end

    for i = env3c.Nz - I:-1:zSharp
        % ３底面からくびれまで
        count = count + 1;
        r_CCS(count) = r(rSharp);
        z_CCS(count) = z(i);
        psi_CCS(count) = psi(rSharp, i);
        Bz_CCS(count) = Bz(rSharp, i);
        Br_CCS(count) = Br(rSharp, i);
    end

    for i = rSharp + 1:env3c.Nr - I - 1
        % くびれからくびれの外側まで
        count = count + 1;
        r_CCS(count) = r(i);
        z_CCS(count) = z(zSharp);
        psi_CCS(count) = psi(i, zSharp);
        Bz_CCS(count) = Bz(i, zSharp);
        Br_CCS(count) = Br(i, zSharp);
    end

    for i = zSharp:-1:zCenter
        % くびれの外側から中央へ
        count = count + 1;
        r_CCS(count) = r(env3c.Nr - I);
        z_CCS(count) = z(i);
        psi_CCS(count) = psi(env3c.Nr - I, i);
        Bz_CCS(count) = Bz(env3c.Nr - I, i);
        Br_CCS(count) = Br(env3c.Nr - I, i);
    end

    if dispFigure

        fh = figure;
        set(fh, 'position', [50 50 600 600], 'Name', ['psi contour'], 'NumberTitle', 'off');

        ah = subplot(3, 1, 1);
        set(ah, 'box', 'on', 'FontSize', fs, 'FontName', FONT);

        hold on
        pcolor(z, r, sqrt(Bz.^2 + Br.^2) * 1e3)
        axis equal
        shading interp

        ca = [-100 100];
        caxis(ca);

        v = linspace(-100, 100, 31);
        contour(z, r, sqrt(Bz.^2 + Br.^2) * 10000, v, 'k')
        v = linspace(0, 5, 11);
        contour(z, r, sqrt(Bz.^2 + Br.^2) * 10000, v, 'k')
        plot([env3c.zmin env3c.zmin env3c.z1 env3c.z2 env3c.z3 env3c.z4 env3c.zmax env3c.zmax], [env3c.rmin env3c.rmirror env3c.rmirror env3c.rmax env3c.rmax env3c.rmirror env3c.rmirror env3c.rmin], 'r', 'LineWidth', 2);
        title('|Bp| [mT]');
        set(ah, 'Box', 'on', 'FontSize', fs, 'FontName', FONT, 'FontWeight', 'bold'); %,'xlim', [-1.2 1.2], 'ylim', [0 0.8]);

        ah = subplot(3, 1, 2);
        set(ah, 'box', 'on', 'FontSize', fs, 'FontName', FONT);

        hold on
        pcolor(z, r, psi * -1000)
        axis equal
        shading interp

        ca = [-5 5];
        caxis(ca);

        v = linspace(-15, 15, 11);
        contour(z, r, psi * -1000, v, 'k')
        v = linspace(-5, 5, 11);
        contour(z, r, psi * -1000, v, 'w')
        %contour(z,r,psi*-1000,-2:0.1:0,'r')
        plot([env3c.zmin env3c.zmin env3c.z1 env3c.z2 env3c.z3 env3c.z4 env3c.zmax env3c.zmax], [env3c.rmin env3c.rmirror env3c.rmirror env3c.rmax env3c.rmax env3c.rmirror env3c.rmirror env3c.rmin], 'r', 'LineWidth', 2);

        % plot(z_PS,r_PS,'r','LineWidth',2);
        % plot(z_LCFS,r_LCFS,'g','LineWidth',2);
        title('poloidal flux [mWb]');
        set(ah, 'Box', 'on', 'FontSize', fs, 'FontName', FONT, 'FontWeight', 'bold'); %,'xlim', [-1.2 1.2], 'ylim', [0 0.8]);

        ah = subplot(3, 1, 3);
        set(ah, 'box', 'on', 'FontSize', fs, 'FontName', FONT);
    end

    %% Calc eddy current
    if 1
        psi_eddy = psi - psi_v_3c + psi222 + psi2222;

        delz = z(2) - z(1);
        delr = r(2) - r(1);

        jt1 = (-psi_eddy(601, 1017:1307) / delr / r(602) ...
            - psi_eddy(601, 1017:1307) / 2 / r(602)^2) ...
            / 2 / pi / 4 / pi / 1e-7/1e6;

        jt2 = -squeeze(psi_eddy(600:-1:502, 1306)) / delz ./ r(600:-1:502)' ...
            / 2 / pi / 4 / pi / 1e-7/1e6;

        jt3 = (-psi_eddy(502, 1306:2032) / delr / r(503) ...
            - psi_eddy(502, 1306:2032) / 2 / r(503)^2) ...
            / 2 / pi / 4 / pi / 1e-7/1e6;

        jt4 = -squeeze(psi_eddy(501:-1:2, 2032)) / delz ./ r(501:-1:2)' ...
            / 2 / pi / 4 / pi / 1e-7/1e6;

        jt5 = (-psi_eddy(2, 2031:-1:1017) / delr / r(1) ... % d^2 psi/dr^2
            + psi_eddy(2, 2031:-1:1017) / 2 / r(1)^2) ... % d psi/dr
            / 2 / pi / 4 / pi / 1e-7/1e6;

        L1 = delz * (1:length(jt1));
        L2 = delz * length(jt1) + delr * (1:length(jt2));
        L3 = delz * length(jt1) + delr * length(jt2) + delz * (1:length(jt3));
        L4 = delz * length(jt1) + delr * length(jt2) + delz * length(jt3) + delr * (1:length(jt4));
        L5 = delz * length(jt1) + delr * length(jt2) + delz * length(jt3) + delr * length(jt4) + delz * (1:length(jt5));
        L6 = delz * length(jt1) + delr * length(jt2) + delz * length(jt3) + delr * length(jt4) + delz * length(jt5) + ...
            delz * (1:length(jt5));
        L7 = delz * length(jt1) + delr * length(jt2) + delz * length(jt3) + delr * length(jt4) + delz * length(jt5) + ...
            delz * length(jt5) + delr * (1:length(jt4));
        L8 = delz * length(jt1) + delr * length(jt2) + delz * length(jt3) + delr * length(jt4) + delz * length(jt5) + ...
            delz * length(jt5) + delr * length(jt4) + delz * (1:length(jt3));
        L9 = delz * length(jt1) + delr * length(jt2) + delz * length(jt3) + delr * length(jt4) + delz * length(jt5) + ...
            delz * length(jt5) + delr * length(jt4) + delz * length(jt3) + delr * (1:length(jt2));
        L10 = delz * length(jt1) + delr * length(jt2) + delz * length(jt3) + delr * length(jt4) + delz * length(jt5) + ...
            delz * length(jt5) + delr * length(jt4) + delz * length(jt3) + delr * length(jt2) + delz * (1:length(jt1));

        Length = [L1 L2 L3 L4 L5 L6 L7 L8 L9 L10];
        jeddy = [jt1 jt2' jt3 jt4' jt5 jt5(end:-1:1) jt4(end:-1:1)' jt3(end:-1:1) jt2(end:-1:1)' jt1(end:-1:1)];
    end

    if dispFigure
        fh = figure;
        hold on
        plot(L1, jt1, 'r');
        plot(L2, jt2, 'r');
        plot(L3, jt3, 'r');
        plot(L4, jt4, 'r');
        plot(L5, jt5, 'r');
        plot(L6, jt5(end:-1:1), 'r');
        plot(L7, jt4(end:-1:1), 'r');
        plot(L8, jt3(end:-1:1), 'r');
        plot(L9, jt2(end:-1:1), 'r');
        plot(L10, jt1(end:-1:1), 'r');
        plot(Length, cumtrapz(Length, jeddy))
        xlabel('Distance (m)');
        ylabel('Eddy current density (MA/m)');
    end

    if saveflag == "origin"
        fp = fopen(save_dir + '/jeddy.txt', 'w');

        for j = 1:length(Length)
            fprintf(fp, '%f\t%f\n', Length(j), jeddy(j));
        end

        fclose(fp);
    elseif saveflag == "modan"
        % jeddy = load("CCS_dataset_gs/UTST_numel_2033/merged").jeddy;
        r_CCS = [r_CCS fliplr(r_CCS)];
        z_CCS = [z_CCS -fliplr(z_CCS)];
        Br_CCS = [Br_CCS -fliplr(Br_CCS)];
        Bz_CCS = [Bz_CCS fliplr(Bz_CCS)];
        psi_CCS = [psi_CCS fliplr(psi_CCS)];
        % save([save_dir + '/merged.mat'], "r_CCS", "z_CCS", "psi_CCS", "Bz_CCS", "Br_CCS", "Length", "jeddy", "psi", "Bz", "Br", "z", "r", "env3c");
        save([save_dir + '/merged.mat'], "r_CCS", "z_CCS", "psi_CCS", "Bz_CCS", "Br_CCS", "psi", "Bz", "Br", "z", "r", "env3c");
    end
end
