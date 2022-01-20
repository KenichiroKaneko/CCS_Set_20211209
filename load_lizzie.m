function [FL, BZ, ExtCOIL, Ip] = load_lizzie(PARAM, CONFIG)
    % バイナリデータをセンサーごとに読み込む
    % ここでは再構成する時間に合わせてデータを切り出さない

    data = UTSTdataload(PARAM.input_file_directory, PARAM.date, PARAM.shotnum, 'NJ');

    if CONFIG.CalibPF3
        [calib_FL, calib_BZ, FL_PF3_t, BZ_PF3_t] = calib_TF3(PARAM);
    end
    if CONFIG.CalibTF
        data_TF_shot = UTSTdataload(PARAM.input_file_directory, PARAM.date, PARAM.shotnum_TF, 'NJ'); 
    end
    t = 0.5:0.5:15000;

    % FLの読み込み
    FL_index = [19:42, 18, 44:53];
    % minus_index = [19, 20, 21, 22, 27, 28, 29, 31, 32, 33, 37, 38, 41, 42, 44, 45, 48, 50, 51, 53];
    minus_index = [5,7,8,12,16,17,18,21,29,34,35];
    FL_num = 35;

    for i = 1:FL_num
        if CONFIG.CalibTF
            FL(i, :) = data(FL_index(i)).sig - data_TF_shot(FL_index(i)).sig;
        else
            FL(i, :) = data(FL_index(i)).sig;
        end
        FL(i, :) = FL(i, :) - mean(FL(i, 1000:2000));
        if sum(i == minus_index)
            FL(i, :) = (-1) * FL(i, :);
        end
    end

    % BZの読み込み
    BZ_index = [54:93];
    % minus_index = [74, 79, 80, 81, 82, 84, 85, 86, 92];
    minus_index = [21,26,27,28,29,31,32,33,39];
    BZ_num = 40;
    NS = 242; % コイルの巻き数x断面積
    freq = fftshift((0:BZ_num-1) * 2000000 / BZ_num) - 1000000;
    lowpass_freq = 1000 * 1000;

    for i = 1:BZ_num
        MP_sig_fft = fft(data(BZ_index(i)).sig);
        MP_sig_fft(freq > lowpass_freq) = 0;
        MP_sig_fft(freq < -lowpass_freq) = 0;
        MP_sig = ifft(MP_sig_fft);
        BZ(i, :) = cumtrapz(MP_sig - mean(MP_sig(1000:2000))) * 0.5e-6 * NS;
        if CONFIG.CalibTF
            MP_sig_TF_fft = fft(data_TF_shot(BZ_index(i)).sig);
            MP_sig_TF_fft(freq > lowpass_freq) = 0;
            MP_sig_TF_fft(freq < -lowpass_freq) = 0;
            MP_sig_TF = ifft(MP_sig_TF_fft);
            BZ_TF(i, :) = cumtrapz(MP_sig_TF - mean(MP_sig_TF(1000:2000))) * 0.5e-6 * NS;
            BZ(i, :) = BZ(i, :) - BZ_TF(i, :);        
        end
        if sum(i == minus_index)
            BZ(i,:) = (-1) * BZ(i, :);
        end
    end

    if CONFIG.CalibPF3
        BZ = BZ .* calib_BZ';
        FL = FL .* calib_FL';
    end
    
    % Ipの読み込み
    Ip = data(2).sig;
    Ip = Ip - mean(Ip(1000:2000));

    % コイル電流の読み込み
    EF_voltage = 120;
    ExtCOIL.NUM = 10;
    ExtCOIL.NAME = ["EFL", "EFU", "PF1L", "PF1U", "PF2L", "PF2U", "PF3L", "PF3U", "PF4L", "PF4U"];
    ExtCOIL.R = [0.80, 0.80, 0.20, 0.20, 0.665, 0.665, 0.750, 0.750, 0.685, 0.685];
    ExtCOIL.Z = [-1.07, 1.07, -1.10, 1.10, -0.80, 0.80, -0.675, 0.675, -0.50, 0.50];
    ExtCOIL.C = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    ExtCOIL.N = [200, 200, 8, 8, 3, 3, 8, 8, 3, 3];
    ExtCOIL.I = [0.28, 0.28, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    ExtCOIL.I(1) =- (0.849 * (1.19 * EF_voltage - 5.32) - 5.56) / 1000;
    ExtCOIL.I(2) =- (0.849 * (1.19 * EF_voltage - 5.32) - 5.56) / 1000;

    % データ（NJ）はUpper→Lowerの順だが、CCSコードの内部ではLower→Upperの順で読み込んでいるので
    % インデックスに気をつけて読み込んでいる
    PF_index = [3, 2, 5, 4, 7, 6, 16, 8];

    for i = 1:8
        ExtCOIL.I_sig(i, :) = smoothdata(data(PF_index(i)).sig, 'gaussian', 200);
    end

    % この辺で信号処理する予定
    %

end

function [calib_FL, calib_BZ, FL_PF3_t, BZ_PF3_t] = calib_TF3(PARAM)
    % それぞれのセンサー位置に作る磁場磁束を再現し、較正に使う
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

    % PF3shotを利用した較正
    data_PF3 = UTSTdataload(PARAM.input_file_directory, PARAM.date_calib_PF3, PARAM.shotnum_PF3, 'NJ');
    t = 0.5:0.5:15000;

    PF3_Z = [-0.675, 0.675];
    PF3_R = [0.75, 0.75];
    PF3_N = [8, 8];
    PF3_index = [7, 6];
    for i = 1:2
        PF3_I(i, :) = smoothdata(data_PF3(PF3_index(i)).sig, 'gaussian', 200);
        % PF3_I(i, :) = data_PF3(PF3_index(i)).sig;
    end

    ECI = PF3_I.* PF3_N' .* 1000;
    RMYU0 = 4 * pi * 1e-7;

    % FLセンサーに作る磁束を計算する
    [PHI, PHIRs, PHIZs, GSTAR, HSTAR, D] = cal_starb(FL_R, FL_Z, PF3_R, PF3_Z, 0, 0, 0);
    psi = (-1) * RMYU0 .* PHI' * ECI;
    FL_index = [19:42, 18, 44:53];
    minus_index = [5,7,8,12,16,17,18,21,29,34,35];
    for i = 1:FL_num
        FL_PF3(i, :) = data_PF3(FL_index(i)).sig;
        FL_PF3(i, :) = FL_PF3(i, :) - mean(FL_PF3(i, 1000:2000));
        if sum(i == minus_index)
            FL_PF3(i, :) = (-1) * FL_PF3(i, :);
        end
    end
    ind_pf3 = [[1,19]; [2,18]; [3,17]; [4,16]; [5,15]; [6,14]; [7,13]; [8,12]; [9,11]; [20,35]; [21,34]; [22,33]; [23,32]; [24,31]; [25,30]; [26,29]; [27,28]];
    figure()
    for i = 1:FL_num
        calib_FL(i) = max(abs(psi(i, 5000:end)))/max(abs(FL_PF3(i, 5000:end)));
    end
    for i = 1:length(ind_pf3)
        subplot(6, 3, i)
        hold on
        % plot(t(5000:end), calib_FL(ind_pf3(i, 1)) * FL_PF3(ind_pf3(i,1), 5000:end))
        % plot(t(5000:end), calib_FL(ind_pf3(i, 2)) * FL_PF3(ind_pf3(i,2), 5000:end))
        plot(t(5000:end), FL_PF3(ind_pf3(i,1), 5000:end))
        plot(t(5000:end), FL_PF3(ind_pf3(i,2), 5000:end))
        plot(t(5000:end), psi(ind_pf3(i, 1), 5000:end))
        plot(t(5000:end), psi(ind_pf3(i, 2), 5000:end))
        legend(num2str(ind_pf3(i,1)),num2str(ind_pf3(i,2)),['cal' num2str(ind_pf3(i,1))],['cal' num2str(ind_pf3(i,2))], 'Location', 'eastoutside')
    end
    
    % BZセンサーに作る磁場を計算する
    [PHI, PHIRs, PHIZs, GSTAR, HSTAR, D] = cal_starb(BZ_R, BZ_Z, PF3_R, PF3_Z, 0, 0, 0);
    BZ = RMYU0 .* PHIRs' * ECI ./ BZ_R';
    BZ_index = [54:93];
    NS = 242;
    minus_index = [21,26,27,28,29,31,32,33,39];
    for i = 1:BZ_num
        BZ_PF3(i, :) = data_PF3(BZ_index(i)).sig;
        BZ_PF3(i, :) = cumtrapz(BZ_PF3(i, :) - mean(BZ_PF3(i, 1000:2000))) * 0.5e-6 * NS;
        if sum(i == minus_index)
            BZ_PF3(i, :) = (-1) * BZ_PF3(i, :);
        end
    end

    ind_pf3 = [[1,18]; [2,17]; [3,16]; [4,15]; [5,14]; [6,13]; [7,12]; [8,11]; [9,10]; [19, 40]; [20,39]; [21,38]; [22,37]; [23,36]; [24,35]; [25,34]; [26,33]; [27,32]; [28,31]; [29,30]];
    figure()
    for i = 1:BZ_num
        calib_BZ(i) = max(abs(BZ(i, 5000:end)))/max(abs(BZ_PF3(i, 5000:end)));
    end
    for i = 1:length(ind_pf3)
        subplot(7,3,i)
        hold on
        % plot(t(5000:end), calib_BZ(ind_pf3(i, 1)) * BZ_PF3(ind_pf3(i,1), 5000:end))
        % plot(t(5000:end), calib_BZ(ind_pf3(i, 2)) * BZ_PF3(ind_pf3(i,2), 5000:end))
        plot(t(5000:end), BZ_PF3(ind_pf3(i,1), 5000:end))
        plot(t(5000:end), BZ_PF3(ind_pf3(i,2), 5000:end))
        plot(t(5000:end), BZ(ind_pf3(i, 1), 5000:end))
        plot(t(5000:end), BZ(ind_pf3(i, 2), 5000:end))
        legend(num2str(ind_pf3(i,1)),num2str(ind_pf3(i,2)),['cal' num2str(ind_pf3(i,1))],['cal' num2str(ind_pf3(i,2))], 'Location', 'eastoutside')
    end

    ccstime = round((str2double(PARAM.time_CCS)) / 0.5);
    BZ_PF3_t = BZ(:, ccstime);
    FL_PF3_t = psi(:, ccstime);
end

function data = UTSTdataload(dir, date, shotnum, type)
    p = join([dir '/utst/' date '/' type shotnum '.hdr'], "");
    p = "./input/utst/210906/NJ007.hdr"; % FLの増幅率をchごとに書いたheaderファイル

    if (exist(p, 'file') == 2)
        file = fopen(p, 'r');
    else
        disp(p);
        data = [];
        return;
    end

    version = fscanf(file, '%d', [1, 1]); % % 1
    num = fscanf(file, '%d', [1, 1]); % % ninum
    data_s = struct('ch', cell(1), 'name', cell(1), 'length', cell(1), 'freq', cell(1), 'offset', cell(1), 'range', cell(1), 'factor', cell(1), ...
        'raw', cell(1), 'sig', cell(1), 'time', cell(1));
    data = repmat(data_s, num, 1);

    for i = 1:num
        data(i).ch = fscanf(file, '%d', [1, 1]); % % ch_number
        data(i).name = fscanf(file, '%s', [1, 1]); % % ch_name
        data(i).length = fscanf(file, '%d', [1, 1]); % % ch_length
        data(i).freq = fscanf(file, '%f', [1, 1]) * 1e-6; % % ch_freq_in_microsec
        data(i).offset = fscanf(file, '%f', [1, 1]); % % ch_offset
        data(i).range = fscanf(file, '%f', [1, 1]); % % ch_fullrange
        data(i).factor = fscanf(file, '%f', [1, 1]); % % ch_factor_to_be_multiplied
    end

    fclose(file);

    p = join([dir '/utst/' date '/' type shotnum '.dat'], "");
    disp(p)

    switch version
        case 1

            if (exist(p, 'file') == 2)
                file = fopen(p, 'r');
            else
                disp(p);
                data = [];
                return;
            end

        case 2

            if (exist(p, 'file') == 2)
                file = fopen(p, 'r', 'b');
            else
                disp('No such a file or directory')
                disp(p)
                data = [];
                return;
            end

    end

    for i = 1:num
        data(i).raw = fread(file, [data(i).length, 1], 'int16');
        data(i).sig = data(i).raw / 2^15 * data(i).range * data(i).factor;
        data(i).time = 0:data(i).freq:data(i).freq * (data(i).length - 1);
    end

    fclose(file);
end

function [Bz, Psi] = EF_calc_for_CCS_probe(r, z, r_size, z_size, V)
    R_EF = 855 * 1e-3;
    z_EF = 1.05;
    Turn_EF = 200;
    mu0 = 4 * pi * 1e-7;
    V = 120;
    I =- (0.849 * (1.19 * V - 5.32) - 5.56);
    Psi = zeros(z_size, r_size);

    for i = 1:1
        alpha_u = sqrt((R_EF + r).^2 + (z - z_EF).^2);
        alpha_l = sqrt((R_EF + r).^2 + (z + z_EF).^2);
        k2_u = (4 * R_EF * r) ./ alpha_u.^2;
        k_u = sqrt(k2_u);
        k2_l = (4 * R_EF * r) ./ alpha_l.^2;
        k_l = sqrt(k2_l);
        m_u = k2_u;
        m_l = k2_l;
        [K_u, E_u] = ellipke(m_u);
        [K_l, E_l] = ellipke(m_l);
        Bz_upper = ((mu0 * Turn_EF * I) ./ (2 * pi)) .* (1 ./ alpha_u) .* (K_u + ((R_EF^2 - r.^2 - (z - z_EF).^2) ./ ((R_EF - r).^2 + (z - z_EF).^2)) .* E_u);
        Bz_lower = ((mu0 * Turn_EF * I) ./ (2 * pi)) .* (1 ./ alpha_l) .* (K_l + ((R_EF^2 - r.^2 - (z + z_EF).^2) ./ ((R_EF - r).^2 + (z + z_EF).^2)) .* E_l);
        Bz = Bz_upper + Bz_lower;

        A_phai_u = (mu0 * Turn_EF * I) ./ (pi) .* sqrt(R_EF ./ r) .* (1 ./ k_u) .* ((1 - k_u.^2 ./ 2) .* K_u - E_u);
        A_phai_l = (mu0 * Turn_EF * I) ./ (pi) .* sqrt(R_EF ./ r) .* (1 ./ k_l) .* ((1 - k_l.^2 ./ 2) .* K_l - E_l);

        psi_u = 2 * pi * r .* A_phai_u;
        psi_l = 2 * pi * r .* A_phai_l;
        Psi = psi_u + psi_l;

    end

end
