% 数値解の作成と、CCS用の入力ファイルの作成

close all; clear all;
% 諸変数の定義
EPS = 1.0e-8; % 許容誤差
Mu = 4 * pi * 1.0e-7; % 真空の透磁率
kappa = 1.38e-23; % ボルツマン定数

% Sw_initを受け取る
param.Sw_init = 1;

% ファイル名や定数の定義
% filenames = ["z201_r101Ip70"];
% filenames = ["z1100_r602", "z1200_r602", "z1300_r602"];
% filenames = ["z1000_r602Ip50", "z1200_r602Ip50"];
% filenames = ["0600Ip50","0640Ip50","0680Ip50","0720Ip50","0760Ip50","0840Ip50",...
%              "0880Ip50","0960Ip50"];
% Nz_list = [600, 640, 680, 720, 760, 840, 880, 960];
% zmax_list = [-0.4088, -0.3695, -0.3302, -0.2909, -0.2516, -0.1730, -0.1337, -0.0550];
% z34_list = [0, 0, 0, 0, 0, 0, 0, 0];

filenames = ["0920Ip50"];
Nz_list = [920];
zmax_list = [-0.0943];
z34_list = [0];

for fileNum = 1:length(filenames)

    % ファイルから値の読み込み
    param.Name = filenames(fileNum);
    param.Nz = 2033;
    param.Nr = 602;
    param.Nz = Nz_list(fileNum);
    param.Nz = 401;
    param.Nr = 201;
    param.Nz_orig = param.Nz;
    param.Nr_orig = param.Nr;
    param.zmin = -0.9985;
    param.zmax = 0.9985; % 2033
    % param.zmax = -0.0157; % 1000
    % param.zmax = -0.0943; % 0900
    % param.zmax = -0.2123; % 0800
    % param.zmax = -0.4088; % 0600
    % param.zmax = zmax_list(fileNum); % 1000
    param.rmin = 0.10815;
    param.rmax = 0.694;
    param.z1 = -0.285;
    param.z2 = -0.285;
    % param.z3 = 0;
    % param.z4 = 0;
    param.z3 = 0.285;
    param.z4 = 0.285;
    param.rmirror = 0.5985;
    param.routine_num = 0;
    param.zlimiter = 0;
    param.rlimiter = 0;
    param.rlimiter = 0.1656627;

    param.err = [];

    % GSplot用
    env3c.Nz = 2033;
    env3c.Nr = 602;
    % env3c.Nz = 201;
    % env3c.Nr = 101;
    env3c.zmin = -0.9985;
    env3c.zmax = 0.9985;
    env3c.rmin = 0.10815;
    env3c.rmax = 0.694;
    env3c.z1 = -0.285;
    env3c.z2 = -0.285;
    env3c.z3 = 0.285;
    env3c.z4 = 0.285;
    env3c.rmirror = 0.5985;
    env3c.delr = (env3c.rmax - env3c.rmin) / (env3c.Nr - 1);
    env3c.delz = (env3c.zmax - env3c.zmin) / (env3c.Nz - 1);
    env3c = InitParams_slow(env3c);
    env3c.psi0 = zeros(env3c.Nz, env3c.Nr);
    [env3c, env3c.psi0] = cal_psi(env3c, env3c.psi0);

    % ファイルから値の読み込み
    param.alpha = 1.01;
    param.beta = 1.30;
    param.P0 = 500; % 0~500~
    param.Iplasma = -50e3; % 大きさ
    param.Itfc = 400e3; % あまり変わらない
    param.errormax = 1e-8;
    param.itermax = 1000;
    param.RBt = Mu * param.Itfc / 2 / pi;

    % 初期状態を設定
    param = define_params(param, param.Nz, param.Nr);

    %% DEFINE OMEGA
    param.omega = 1.945;

    %% 荒いメッシュで計算
    [param, psi, mu_jt, p, Ip] = sub_routine_init(param, 4, 4);
    [param, psi, mu_jt, p, Ip] = sub_routine(param, 4, 4, psi, mu_jt, p, Ip);
    % [param, psi, mu_jt, p, Ip] = sub_routine_init(param, 2, 2);
    [param, psi, mu_jt, p, Ip] = sub_routine(param, 2, 2, psi, mu_jt, p, Ip);
    % [param, psi, mu_jt, p, Ip] = sub_routine_init(param, 1.5, 1.5);
    [param, psi, mu_jt, p, Ip] = sub_routine(param, 1.5, 1.5, psi, mu_jt, p, Ip);
    % %% メインのルーチン用に変形
    param.Nz = param.Nz_orig;
    param.Nr = param.Nr_orig;
    psi = imresize(psi, [param.Nz, param.Nr]);
    % psi = zeros(param.Nz, param.Nr);
    mu_jt = zeros(param.Nz, param.Nr);
    p = zeros(param.Nz, param.Nr);
    Ip = zeros(param.Nz, param.Nr);
    param = define_params(param, param.Nz, param.Nr);
    % [param, psi] = cal_initial_flux(param, psi);
    psi = psi_round_revive(psi, param);

    %% メインのルーチン
    for i = 1:1000000
        [psi_bar, mu_jt, p, Ip] = cal_jt(param, psi, mu_jt, p, Ip);
        format long
        [error, psi] = cal_flux_C(param, psi, mu_jt);


        if (rem(i, 100) == 0)
            disp([num2str(i) ':' num2str(error)])
        end

        if error < param.errormax
            disp(['Convergence! Iteration number is ' num2str(param.routine_num)]);
            break
        else
            param.routine_num = param.routine_num + 1;
            param.err(param.routine_num) = error;
        end

    end

    psi0 = zeros(param.Nz, param.Nr);
    [param, psi0] = cal_psi(param, psi0);

    % 変数の保存
    % save(output_dir + "vars", "param", "psi", "mu_jt", "Ip", "p", "psi0");
    vars.psi = psi;
    vars.mu_jt = mu_jt;
    vars.Ip = Ip;
    vars.p = p;
    vars.psi0 = psi0;

    %% グラフ
    r = param.rmin + param.delr * (0:param.Nr - 1);
    z = param.zmin + param.delz * (0:param.Nz - 1);
    [rr zz] = meshgrid(r, z);
    v = linspace(-20, 20, 21);
    figure()
    hold on
    contour(rr, zz, psi * 1000, v, 'r')
    param.r_limiterpos;
    param.z_limiterpos;
    plot(r(param.r_limiterpos+1), z(param.z_limiterpos+1), 'o');
    axis equal

    figure()
    semilogy(param.err)

    % save('vars_gs_temp');
    % gs2ccs(param, env3c, vars);

end


%% psiを更新するルーチン。
function [error, psi] = cal_flux_C(param, psi, mu_jt)
    error = 0;
    vec3 = 1 / (param.delz^2);
    vec4 = 1 / (param.delz^2);

    for j = 2:param.Nr - 1
        radius = param.rmin + (j - 1) * param.delr;
        vec1 = 1 / (param.delr^2) * radius / (radius + 0.5 * param.delr);
        vec2 = 1 / (param.delr^2) * radius / (radius - 0.5 * param.delr);
        vec5 = vec1 + vec2 + vec3 + vec4;

        for i = 2:param.Nz - 1

            if param.psi_inside(i, j) == 1
                psi_old = psi(i, j);
                psi(i, j) = param.omega * ( ...
                    vec1 * psi(i, j + 1) + ...
                    vec2 * psi(i, j - 1) + ...
                    vec3 * psi(i + 1, j) + ...
                    vec4 * psi(i - 1, j) + ...
                    2 * pi * radius * mu_jt(i, j)) / vec5 + ...
                    (1 - param.omega) * psi_old;

                if (abs(psi(i, j) - psi_old) > error)
                    error = abs(psi(i, j) - psi_old);
                end

            end

        end

    end

end

%% ψの初期状態を置く
function [param, psi] = cal_psi(param, psi)
    Mu = 4 * pi * 1.0e-7; % 真空の透磁率
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

end

%% psiの初期値を適当に決める
function [param, psi] = cal_initial_condition(param, psi)
    psi_axis = -1.75;
    psi_edge = 0;
    zc = (param.Nz - 1) / 2 * param.delz;
    rc = (param.Nr - 1) / 3 * param.delr;
    mtrxz = (0:1:param.Nz - 1)' .* param.delz;
    mtrxr = (0:1:param.Nr - 1) .* param.delr;
    radius = sqrt(((mtrxz - zc) ./ zc).^2 * ones(1, param.Nr) + ones(param.Nz, 1) * (mtrxr ./ rc).^2);
    psi_b = (psi_axis - psi_edge) * (1 - radius);
    % 2020/12/03
    psi_b(psi_b > 0) = 0;
    psi = psi - psi .* param.psi_inside + psi_b .* param.psi_inside;
end

%% 磁束の値から電流jtを求める
function [psi_bar, mu_jt, p, Ip] = cal_jt(param, psi, mu_jt, p, Ip)

    EPS = 1.0e-8; % 許容誤差
    Mu = 4 * pi * 1.0e-7; % 真空の透磁率

    % psi_axis = min(psi(param.psi_inside == 1));
    psi_axis = min(min(psi(param.psi_inside == 1)));
    psi_edge = min(psi(param.psi_round == 1));

    if (param.z_limiterpos ~= 1 || param.r_limiterpos ~= 1)
        if (psi(param.z_limiterpos, param.r_limiterpos) < psi_edge)
            psi_edge = psi(param.z_limiterpos, param.r_limiterpos);
        end
    end

    % psi_edge = min(psi(66, :));

    if (psi_edge <= (psi_axis) + EPS)
        error("psi_edge <= psi_axis ...");
    end

    radius = param.rmin + ones(param.Nz, 1) * (0:1:param.Nr - 1) .* param.delr;
    psi_resion = param.psi_inside;
    psi_resion(psi > psi_edge | psi < psi_axis) = 0;
    psi_deff = psi_axis - psi_edge;
    ip = CALGP_t(param, (psi - psi_edge)/ psi_deff) / psi_deff;

    aaa = sum(((psi - psi_edge) .* psi_resion ./ psi_deff).^(2 * param.beta - 1) ./ radius ./ psi_deff, "all");
    bbb = sum(((psi - psi_edge) .* psi_resion ./ psi_deff).^(param.beta - 1) ./ radius ./ psi_deff, "all");
    ccc = sum(radius .* ip .* psi_resion, "all");

    aaa = aaa * (2 * pi * param.RBt^2 * param.beta * param.delr * param.delz);
    bbb = bbb * (-2 * pi * param.RBt^2 * param.beta * param.delr * param.delz);
    ccc = ccc * (2 * pi * Mu * param.delr * param.delz);
    ccc = ccc - param.Iplasma * Mu;

    if abs(aaa) < EPS
        param.gamma = -ccc / bbb;
    else
        abc = [aaa, bbb, ccc];
        sol = min(roots(abc));

        if isreal(sol)
            param.gamma = sol;
        else
            figure()
            contour(psi)
            error("Solutions are complex.");
        end

    end

    radius = param.rmin + ones(param.Nz, 1) * (0:1:param.Nr - 1) .* param.delr;

    mu_jt = 2 * pi .* (radius .* Mu ...
        .* CALGP_t(param, (psi - psi_edge) ./ psi_deff) ./ psi_deff ...
        + CALIp_t(param, (psi - psi_edge) ./ psi_deff) ...
        .* CALGIp_t(param, (psi - psi_edge) ./ psi_deff) ./ psi_deff ./ radius);
    p = CALp_t(param, (psi - psi_edge) ./ psi_deff);
    Ip = CALIp_t(param, (psi - psi_edge) ./ psi_deff);
    mu_jt = mu_jt .* param.psi_inside;
    p = p .* param.psi_inside;
    Ip = Ip .* param.psi_inside;
    psi_bar = psi_axis / psi_edge;
end

function pp = CALp_t(param, psibar)
    pp = zeros(param.Nz, param.Nr);
    pp(psibar >= 0) = param.P0 .* power(psibar(psibar >= 0), param.alpha);
    pp(psibar < 0) = 0;
end

function gpp = CALGP_t(param, psibar)
    gpp = zeros(param.Nz, param.Nr);
    gpp(psibar >= 0) = param.alpha .* param.P0 .* power(psibar(psibar >= 0), param.alpha - 1);
    gpp(psibar < 0) = 0;
end

function ip = CALIp_t(param, psibar)
    ip = zeros(param.Nz, param.Nr);
    ip(psibar >= 0) = param.RBt .* (1 - param.gamma .* power(psibar(psibar >= 0), param.beta));
    ip(psibar < 0) = param.RBt;
end

function gip = CALGIp_t(param, psibar)
    gip = zeros(param.Nz, param.Nr);
    gip(psibar >= 0) = (-param.RBt * param.gamma * param.beta) .* power(psibar(psibar >= 0), param.beta - 1);
    gip(psibar < 0) = 0;
end

%% psiを更新するルーチン,配列を利用したものだがうまく動かない
function [error, psi] = cal_flux(param, psi, mu_jt)
    vec3 = 1 / (param.delz^2);
    vec4 = 1 / (param.delz^2);

    radius = param.rmin + ones(param.Nz, 1) * (0:1:param.Nr - 1) .* param.delr;
    vec1 = 1 / (param.delr^2) .* radius ./ (radius + 0.5 * param.delr);
    vec2 = 1 / (param.delr^2) .* radius ./ (radius - 0.5 * param.delr);

    vec5 = vec1 + vec2 + vec3 + vec4;

    psi_old = psi;
    psi_new = param.omega * ( ...
        vec1 .* circshift(psi, [0 -1]) ...
        +vec2 .* circshift(psi, [0 1]) ...
        +vec3 .* circshift(psi, [1 0]) ...
        +vec4 .* circshift(psi, [-1 0]) ...
        +2 * pi .* radius .* mu_jt) ./ vec5 ...
        + (1 - param.omega) * psi_old;

    psi = psi_old .* param.psi_round + psi_new .* param.psi_inside;

    C = max(abs(psi(param.psi_inside == 1) - psi_old(param.psi_inside == 1)));
    error = C;
end

%%
function mu_jt = addPFcurrent(param, psi, mu_jt)
    Mu = 4 * pi * 1.0e-7; % 真空の透磁率
    j = fix((param.coil_r(1, :) - param.rmin) ./ param.delr);
    r1 = param.coil_r(1, :) - (param.rmin + j * param.delr);
    i = fix((param.coil_z(1, :) - param.zmin) ./ param.delz);
    z1 = param.coil_z(1, :) - (param.zmin + i * param.delz);
    cur = Mu * param.coil_Ic(:)' ./ (param.delr^2 * param.delz^2);

    mu_jt(indexToMatrix(param, i + param.Nz .* j)) = mu_jt(indexToMatrix(param, i + param.Nz .* j)) + cur .* (param.delr - r1) .* (param.delz - z1);
    mu_jt(indexToMatrix(param, i + param.Nz .* (j + 1))) = mu_jt(indexToMatrix(param, i + param.Nz .* (j + 1))) + cur .* (r1) .* (param.delz - z1);
    mu_jt(indexToMatrix(param, i + 1 + param.Nz .* j)) = mu_jt(indexToMatrix(param, i + 1 + param.Nz .* j)) + cur .* (param.delr - r1) .* (z1);
    mu_jt(indexToMatrix(param, i + 1 + param.Nz .* (j + 1))) = mu_jt(indexToMatrix(param, i + 1 + param.Nz .* (j + 1))) + cur .* (r1) .* (z1);

    function [z, r] = indexToMatrix(param, num)
        z = fix(num ./ param.Nr);
        r = rem(num, param.Nr);

        if r == 0
            r = param.Nr;
        end

    end

end

% コマンドから実行する関数を選択するための関数
function [param, psi] = cal_initial_flux(param, psi)

    switch param.Sw_init
        case 1
            param = InitParams_slow(param);
            [param, psi] = cal_psi(param, psi);
            disp("after cal_psi");
            [param, psi] = cal_initial_condition(param, psi);
            disp("after cal_initial_condition")
        % case 2
        %     ReadIn_data(psi, output_dir + "flux0.txt");
        %     cal_initial_condition(param, psi);
        % case 3
        %     ReadIn_data(psi, output_dir + "flux.txt");
        otherwise
            disp("Input correct number.");
    end

end

% コイルに関するデータの読み込み
function param = InitParams_slow(param)
    param.coil_turn = [200 200 8 8 3 3 8 8 3 3];
    param.coil_Ic = [280 280 0 0 0 0 0 0 0 0];
    param.coil_Ic = param.coil_Ic .* param.coil_turn;
    param.coil_z = [-1.07 1.07 -1.1 1.1 -0.8 -0.8 -0.675 0.675 -0.5 0.5];
    param.coil_r = [0.8 0.8 0.2 0.2 0.685 0.685 0.75 0.75 0.685 0.685];
    param.ncoil = length(param.coil_turn);
end

function ReadIn_data(psi, filename)
    psi = readmatrix(filename);
end

function WriteOut_data(psi, filename)
    writematrix(psi, filename, 'Delimiter', 'tab');
end

function param = define_params(param, Nz, Nr)
    param.Nz = Nz;
    param.Nr = Nr;

    % メッシュ1ステップの長さ
    param.delr = (param.rmax - param.rmin) / (Nr - 1);
    param.delz = (param.zmax - param.zmin) / (Nz - 1);

    % zに応じて格子点がいくつあるかrnumに格納
    param.rnum = zeros(1, Nz);

    for i = 1:Nz

        if ((param.zmin + (i - 1) * param.delz) < param.z1 || param.z4 < (param.zmin + (i - 1) * param.delz))
            param.rnum(1, i) = ceil((Nr - 1) * (param.rmirror - param.rmin) / (param.rmax - param.rmin));
        elseif (param.z2 < (param.zmin + (i - 1) * param.delz) && (param.zmin + (i - 1) * param.delz) < param.z3)
            param.rnum(1, i) = ceil(Nr);
        elseif (param.z1 < (param.zmin + (i - 1) * param.delz) && (param.zmin + (i - 1) * param.delz) < param.z2)
            param.rnum(1, i) = ceil((Nr - 1) * ((param.rmirror - param.rmin) / (param.rmax - param.rmin) + ...
                (param.rmax - param.rmirror) / (param.rmax - param.rmin) * (param.zmin + (i - 1) * param.delz - param.z1) / (param.z2 - param.z1)));
        elseif (param.z3 < (param.zmin + (i - 1) * param.delz) && (param.zmin + (i - 1) * param.delz) < param.z4)
            param.rnum(1, i) = ceil((Nr - 1) * (1 - (param.rmax - param.rmirror) / (param.rmax - param.rmin) * ...
                (param.zmin + (i - 1) * param.delz - param.z3) / (param.z4 - param.z3)));
        end

    end

    %% DEFINE REGIONS
    param.psi_inside = zeros(Nz, Nr);

    for i = 2:Nz - 1
        param.psi_inside(i, 2:param.rnum(1, i) - 1) = 1;
    end

    param.psi_round = zeros(Nz, Nr);

    for i = 1:Nz
        param.psi_round(i, 1:param.rnum(1, i)) = 1;
    end

    param.psi_round = param.psi_round - param.psi_inside;

    for i = 2:Nr

        for j = 1:Nz - 1

            if (param.psi_inside(j, i) == 0 && param.psi_inside(j + 1, i) == 1)
                param.psi_round(j, i) = 1;
            elseif (param.psi_inside(j, i) == 1 && param.psi_inside(j + 1, i) == 0)
                param.psi_round(j + 1, i) = 1;
            end

        end

    end

    param.psi_whole = param.psi_inside + param.psi_round;

    % リミターの位置
    if (param.rlimiter < param.rmin || param.rmax < param.rlimiter || param.zlimiter < param.zmin || param.zmax < param.zlimiter)
        % メッシュ外にリミターがあるなら(0,0)にする
        param.r_limiterpos = 1;
        param.z_limiterpos = 1;
    else
        % メッシュ内にリミターがあるなら最も近い格子点を探し，そこを指定する
        r = param.rmin + param.delr * (0:param.Nr - 1);
        z = param.zmin + param.delz * (0:param.Nz - 1);
        [m param.z_limiterpos] = min(abs(z - param.zlimiter));
        [m param.r_limiterpos] = min(abs(r - param.rlimiter));
        % mtrxz = (0:1:Nz - 1)' .* param.delz + param.zmin;
        % mtrxr = (0:1:Nr - 1) .* param.delr + param.rmin;
        % dot = (mtrxz - param.zlimiter).^2 * ones(1, Nr) + ones(Nz, 1) * (mtrxr - param.rlimiter).^2;
        % % 以下は2次元配列の最小値の座標を求めていて，Cに最小値を代入している
        % [C, I] = min(dot(:));
        % % Cは0オリジンだから値が1ずれる．
        % [param.z_limiterpos, param.r_limiterpos] = ind2sub(size(dot), I);
        [a,b] = min(abs(r - param.rlimiter));
    end

end

function [param, psi, mu_jt, p, Ip] = sub_routine_init(param, divz, divr)
    param.Nz = round(param.Nz_orig / divz);
    param.Nr = round(param.Nr_orig / divr);
    param = define_params(param, param.Nz, param.Nr);
    psi = zeros(param.Nz, param.Nr);
    mu_jt = zeros(param.Nz, param.Nr);
    p = zeros(param.Nz, param.Nr);
    Ip = zeros(param.Nz, param.Nr);
    [param, psi] = cal_initial_flux(param, psi);
end

function [param, psi, mu_jt, p, Ip] = sub_routine(param, divz, divr, psi, mu_jt, p, Ip)
    %% サブのルーチン　収束を早める
    param.Nz = round(param.Nz_orig / divz);
    param.Nr = round(param.Nr_orig / divr);
    param = define_params(param, param.Nz, param.Nr);
    disp(['    sub routine param.Nz:' num2str(param.Nz) ' param.Nr:' num2str(param.Nr)])

    psi = imresize(psi, [param.Nz, param.Nr]);
    psi = psi_round_revive(psi, param);
    mu_jt = imresize(mu_jt, [param.Nz, param.Nr]);
    p = imresize(p, [param.Nz, param.Nr]);
    Ip = imresize(Ip, [param.Nz, param.Nr]);

    for i = 1:100000
        [psi_bar, mu_jt, p, Ip] = cal_jt(param, psi, mu_jt, p, Ip);
        %     mu_jt = addPFcurrent(param,psi,mu_jt);
        format long
        [error, psi] = cal_flux_C(param, psi, mu_jt);

        if (rem(i, 100) == 0)
            disp(['    ' num2str(i) ':' num2str(error)])
        end

        if error < param.errormax
            disp(['  Sub routine finished at ' num2str(param.routine_num)]);
            break
        else
            param.routine_num = param.routine_num + 1;
            param.err(param.routine_num) = error;
        end

    end

end

% psi_roundを復活させる関数
function psi = psi_round_revive(psi, param)
    psi_b = psi;
    psi = zeros(param.Nz, param.Nr);
    [param, psi] = cal_psi(param, psi);
    psi = psi .* param.psi_round + psi_b .* (1 - param.psi_round);
end

function gs2ccs(env, env3c, vars)
    % GSplot_CCS_for_finemesh_merge2('z2033_r602')みたいに実行する
    % CCS03cはUTST装置全体の解
    % CCS05は（途中に仮想壁を置いた）部分的な解
    % 部分的な解を全体に拡張する際に、CCS03cのメッシュデータが必要
    % fluxは平衡解、flux0はpsiの初期状態
    disp("gs2ccs start")
    whos
    
    saveflag = "modan"; % origin or modan or 0
    dispFigure = 0;
    realFlag = 0;
    % 保存するdir
    save_dir = "input3/UTST_numel_" + env.Name;
    if not(exist(save_dir, 'dir'))
        mkdir(save_dir);
    end

    % 全体の解のデータ
    flux3c = env3c.psi0;
    psi_v_3c = flux3c';

    % 容器の大きさに拡張したい部分的な解のデータ
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
    if 0
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
