function [FL, BZ, BR, ExtCOIL, REF, JEDDY, POS] = load_num_sol(PARAM, CONFIG)
    % GSコードで作成した数値解を読み込む
    % GSコードでは容器内の磁束しか求められないので、外側を計算して付け足している

    vars = load([PARAM.input_file_directory PARAM.num_sol_name '/merged.mat']);
    env = vars.env3c;

    % プラズマ無しの容器内の磁束
    % psi602 = zeros(env.Nz, env.Nr_orig);
    psi602 = zeros(2033, 602);
    [env, psi602] = cal_psi(env, psi602);
    psi602 = psi602.';
    % プラズマ無しの容器外も含めた磁束
    env.Nr = 800;
    env.rmax = env.rmin + env.delr * (env.Nr - 1);
    psi800 = zeros(env.Nz, env.Nr);
    [env, psi800] = cal_psi(env, psi800);
    psi800 = psi800.';
    % プラズマ有りの容器内の磁束
    psi_orig = vars.psi;
    % プラズマ有りの容器外も含めた磁束（外側に磁束を付け足し）
    % padding = zeros(env.Nr - env.Nr_orig, env.Nz);
    padding = zeros(env.Nr - 602, env.Nz);
    psi602 = [psi602; padding];
    psi_orig = [psi_orig; padding];
    psi = psi_orig - psi602 + psi800; % 磁束の分布

    % % 外側にセンサーをおかない
    % psi = vars.psi;

    z = linspace(env.zmin, env.zmax, env.Nz);
    r = linspace(env.rmin, env.rmax, env.Nr);

    % 磁場の分布の読み込み
    br = vars.Br;
    bz = vars.Bz;

    % プロット用
    REF.Flux = psi' ./ (2 * pi);
    REF.R = r;
    REF.Z = z;
    % 渦電流分布
    vars = load([PARAM.input_file_directory '/UTST_numel_2033/merged.mat']);
    JEDDY(:, 1) = vars.Length;
    JEDDY(:, 2) = vars.jeddy;

    % 各センサー位置での信号を取得
    senpos_B = fileread(PARAM.temporary_file_directory + PARAM.MP_pos);
    senpos_B = strsplit(senpos_B, {'\n', '\t', '\r'});
    senpos_FL = fileread(PARAM.temporary_file_directory + PARAM.FL_pos);
    senpos_FL = strsplit(senpos_FL, {'\n', '\t', '\r'});
    BZ_num = floor((length(senpos_B) - 5) / 5);
    FL_num = floor((length(senpos_FL) - 5) / 5);
    BZ_R = str2double(senpos_B(6:5:5 * BZ_num + 5));
    BZ_Z = str2double(senpos_B(7:5:5 * BZ_num + 5));
    FL_R = str2double(senpos_FL(6:5:5 * FL_num + 5));
    FL_Z = str2double(senpos_FL(7:5:5 * FL_num + 5));

    if CONFIG.RevSenPos
        % z=0で線対称のセンサー配置にする
        BZ_R = [BZ_R fliplr(BZ_R(2:end-1))];
        BZ_Z = [BZ_Z -fliplr(BZ_Z(2:end-1))];
        FL_R = [FL_R fliplr(FL_R(2:end-1))];
        FL_Z = [FL_Z -fliplr(FL_Z(2:end-1))];
        BZ_num = length(BZ_R);
        FL_num = length(FL_R);
    end

    for i = 1:BZ_num
        [m, I_R] = min(abs(r - BZ_R(i)));
        [m, I_Z] = min(abs(z - BZ_Z(i)));
        BR(i) = br(I_R, I_Z);
        BZ(i) = bz(I_R, I_Z);
        POS.Br(i) = r(I_R); POS.Bz(i) = z(I_Z);
    end

    for i = 1:FL_num
        [m, I_R] = min(abs(r - FL_R(i)));
        [m, I_Z] = min(abs(z - FL_Z(i)));
        FL(i) = psi(I_R, I_Z);
        POS.FLr(i) = r(I_R); POS.FLz(i)= z(I_Z);
    end

    % コイル関係
    ExtCOIL.NUM = 10;
    ExtCOIL.NAME = ["EFL", "EFU", "PF1L", "PF1U", "PF2L", "PF2U", "PF3L", "PF3U", "PF4L", "PF4U"];
    ExtCOIL.R = [0.80, 0.80, 0.20, 0.20, 0.665, 0.665, 0.750, 0.750, 0.685, 0.685];
    ExtCOIL.Z = [-1.07, 1.07, -1.10, 1.10, -0.80, 0.80, -0.675, 0.675, -0.50, 0.50];
    ExtCOIL.C = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    ExtCOIL.N = [200, 200, 8, 8, 3, 3, 8, 8, 3, 3];
    ExtCOIL.I = [0.28, 0.28, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

end

function [param, psi] = cal_psi(param, psi)

    Mu = 4 * pi * 1.0e-7;
    mtrxz = (0:1:param.Nz - 1)' .* param.delz + param.zmin;
    mtrxr = (0:1:param.Nr - 1) .* param.delr + param.rmin;

    for k = 1:param.ncoil
        % (z*rの大きさの配列を作成)
        tmp1 = (param.coil_z(1, k) - mtrxz).^2 * ones(1, param.Nr) + ones(param.Nz, 1) * (param.coil_r(1, k) + mtrxr).^2;
        kk = (4 .* ones(param.Nz, 1) * mtrxr .* param.coil_r(1, k)) ./ tmp1;
        [K, E] = ellipke(kk);
        psi = psi + 2 * pi * Mu * param.coil_Ic(1, k) ./ (pi * (kk).^0.5) .* sqrt(param.coil_r(1, k) .* mtrxr) .* ...
            ((1 - kk / 2) .* K - E);
    end
    % save("cal_psi","psi")

end
