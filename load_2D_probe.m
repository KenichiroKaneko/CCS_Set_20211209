function REF = load_2D_probe(PARAM, ExtCOIL)

    data_NI = UTSTdataload(PARAM.input_file_directory, PARAM.date, PARAM.shotnum, 'NI'); % Mz
    data_ICS = UTSTdataload(PARAM.input_file_directory, PARAM.date, PARAM.shotnum, 'ICS'); % Mt

    EF.coilnum = ExtCOIL.NUM;
    EF.r = ExtCOIL.R;
    EF.z = ExtCOIL.Z;
    EF.N = ExtCOIL.N;
    EF.current = ExtCOIL.I;

    % 磁束面が欲しい時間（ステップ数）
    t = round((str2double(PARAM.time_CCS)) / 0.5);
    t = 20000;


    for i = 1:length(data_NI)
        data_NI(i).time(end - 24:end) = [];
        data_NI(i).raw(end - 24:end) = [];
        data_NI(i).sig(end - 24:end) = [];
        offset = mean(data_NI(i).sig(50:round(450 / data_NI(i).freq / 1e6)));
        data_NI(i).sigint = cumtrapz(data_NI(i).sig - offset) * data_NI(i).freq;
    end

    for i = 1:length(data_ICS)
        %    data_ICS(i).time=data_ICS(i).time-0.5e-4; % 50usec delay
        data_ICS(i).time(1:25) = [];
        data_ICS(i).raw(1:25) = [];
        data_ICS(i).sig(1:25) = [];
        offset = mean(data_ICS(i).sig(50:round(450 / data_ICS(i).freq / 1e6)));
        data_ICS(i).sigint = cumtrapz(data_ICS(i).sig - offset) * data_ICS(i).freq;
    end

    ind_Mz = 1:81;
    ind_Mt = [2:82];
    ind_dead_Mt = [66, 73, 74];

    Mzarray = data_NI(ind_Mz);
    Mtarray = data_ICS(ind_Mt);

    r = [0.11 0.18 0.25 0.32 0.39 0.46 0.53 0.60 0.67 ...
        0.11 0.18 0.25 0.32 0.39 0.46 0.53 0.60 0.67 ...
        0.11 0.18 0.25 0.32 0.39 0.46 0.53 0.60 0.67 ...
        0.11 0.18 0.25 0.32 0.39 0.46 0.53 0.60 0.67 ...
        0.11 0.18 0.25 0.32 0.39 0.46 0.53 0.60 0.67 ...
        0.11 0.18 0.25 0.32 0.39 0.46 0.53 0.60 0.67 ...
        0.11 0.18 0.25 0.32 0.39 0.46 0.53 0.60 0.67 ...
        0.11 0.18 0.25 0.32 0.39 0.46 0.53 0.60 0.67 ...
        0.11 0.18 0.25 0.32 0.39 0.46 0.53 0.60 0.67];
    z = [0.23 0.23 0.23 0.23 0.23 0.23 0.23 0.23 0.23 ...
        0.18 0.18 0.18 0.18 0.18 0.18 0.18 0.18 0.18 ...
        0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 0.12 ...
        0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06 ...
        0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 ...
        -0.06 -0.06 -0.06 -0.06 -0.06 -0.06 -0.06 -0.06 -0.06 ...
        -0.12 -0.12 -0.12 -0.12 -0.12 -0.12 -0.12 -0.12 -0.12 ...
        -0.18 -0.18 -0.18 -0.18 -0.18 -0.18 -0.18 -0.18 -0.18 ...
        -0.23 -0.23 -0.23 -0.23 -0.23 -0.23 -0.23 -0.23 -0.23];

    % 角度校正用のテキストファイル
    angle = load([PARAM.input_file_directory, '/utst/', PARAM.date, '/angle2.txt']);
    calib = load('./temp/RBtcalib.mat')

    for i = 1:81
        magdata(i).Btm = Mtarray(i).sigint;
        magdata(i).Btd = Mtarray(i).sig;
        magdata(i).Bttime = Mtarray(i).time;
        magdata(i).Bzm = Mzarray(i).sigint;
        magdata(i).Bzd = Mzarray(i).sig;
        magdata(i).Bztime = Mzarray(i).time;

        if length(angle) > 81
            magdata(i).angle = angle(i + 64, 2) / 180 * pi;
        else
            magdata(i).angle = angle(i, 2) / 180 * pi; %no NULL probe
        end

        magdata(i).r = r(i);
        magdata(i).z = z(i);
        magdata(i).type = 'MID';
    end

    for i = 1:81
        magdata(i).Bz = cos(magdata(i).angle) * magdata(i).Bzm - sin(magdata(i).angle) * magdata(i).Btm;
        magdata(i).Bt = sin(magdata(i).angle) * magdata(i).Bzm + cos(magdata(i).angle) * magdata(i).Btm;
        magdata(i).Bzdev = cos(magdata(i).angle) * magdata(i).Bzd - sin(magdata(i).angle) * magdata(i).Btd;
        magdata(i).Btdev = sin(magdata(i).angle) * magdata(i).Bzd + cos(magdata(i).angle) * magdata(i).Btd;
        magdata(i).Bt2 = magdata(i).Bt / (magdata(i).Bt(floor(end / 2)) - magdata(i).Bt(end)) * (calib.RBt(floor(end / 2)) - calib.RBt(end)) / magdata(i).r;
        magdata(i).Bt2 = magdata(i).Bt2 - magdata(i).Bt2(end) + calib.RBt(end) / magdata(i).r;
        EFBz = 0;

        for j = 1:EF.coilnum
            k = sqrt(4 * EF.r(j) * r(i) ./ ((EF.r(j) + r(i)).^2 + (EF.z(j) - z(i)).^2));
            [K, E] = ellipke(k.^2);
            % EFBz=EFBz+4*pi*1e-7*sqrt((EF.r(j)+r(i)).^2+(EF.z(j)-z(i)).^2).*((1-k.^2/2).*K-E)*EF.N(j)*EF.current(j);
            EFBz = EFBz + (4 * pi * 1e-7/2 / pi) / sqrt((EF.r(j) + r(i))^2 + (EF.z(j) - z(i))^2) * (K + (EF.r(j)^2 - r(i)^2 - (EF.z(j) - z(i))^2) / ((EF.r(j) - r(i))^2 + (EF.z(j) - z(i))^2) * E) * EF.N(j) * EF.current(j);
        end

        magdata(i).Bz = magdata(i).Bz + EFBz;
    end


    REF = 0;

    for i = 1:81
        Bz(i) = magdata(i).Bz(t);
    end
    whos

    % Z = reshape(z, [9, 9]);
    % R = reshape(r, [9, 9]);
    R = meshgrid(linspace(r(1), r(end), 100));
    Z = meshgrid(linspace(z(1), z(end), 100));
    [R, Z] = ndgrid(linspace(r(1), r(end), 30), linspace(z(1), z(end), 30));
    F = scatteredInterpolant([r; z]', Bz', 'linear');
    Bz = F(R, Z);

    % Bz = reshape(Bz, [9, 9]);
    % R = reshape(r, [9, 9]);
    % Z = reshape(z, [9, 9]);
    integrand = 2*pi*Bz;
    flux = cumtrapz(integrand, 1)*(R(2)-R(1));
    

    figure()
    contour( R, Z, flux)
end


function data = UTSTdataload(dir, date, shotnum, type)
    p = join([dir '/utst/' date '/' type shotnum '.hdr'], "");

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
