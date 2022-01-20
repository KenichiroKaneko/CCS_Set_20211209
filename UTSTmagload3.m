function magdata = UTSTmagload3(dir, date, shotnum, EF, probe)

    EF.coilnum = ExtCOIL.NUM;
    EF.r = ExtCOIL.R;
    EF.z = ExtCOIL.Z;
    EF.N = ExtCOIL.N;
    EF.current = EF.I;


    %% load data
    data_NI = UTSTdataload(dir, date, shotnum, 'NI');

    for i = 1:length(data_NI)
        data_NI(i).time(end - 24:end) = [];
        data_NI(i).raw(end - 24:end) = [];
        data_NI(i).sig(end - 24:end) = [];
        offset = mean(data_NI(i).sig(50:round(450 / data_NI(i).freq / 1e6)));
        data_NI(i).sigint = cumtrapz(data_NI(i).sig - offset) * data_NI(i).freq;
    end

    data_NJ = UTSTdataload(dir, date, shotnum, 'NJ');

    for i = 1:length(data_NJ)
        data_NJ(i).time(end - 24:end) = [];
        data_NJ(i).raw(end - 24:end) = [];
        data_NJ(i).sig(end - 24:end) = [];
        offset = mean(data_NJ(i).sig(50:round(450 / data_NJ(i).freq / 1e6)));
        data_NJ(i).sigint = cumtrapz(data_NJ(i).sig - offset) * data_NJ(i).freq;
    end

    data_NK = UTSTdataload(dir, date, shotnum, 'NK');

    for i = 1:length(data_NK)
        data_NK(i).time(end - 24:end) = [];
        data_NK(i).raw(end - 24:end) = [];
        data_NK(i).sig(end - 24:end) = [];
        offset = mean(data_NK(i).sig(50:round(450 / data_NK(i).freq / 1e6)));
        data_NK(i).sigint = cumtrapz(data_NK(i).sig - offset) * data_NK(i).freq;
    end

    data_ICS = UTSTdataload(dir, date, shotnum, 'ICS');

    for i = 1:length(data_ICS)
        %    data_ICS(i).time=data_ICS(i).time-0.5e-4; % 50usec delay
        data_ICS(i).time(1:25) = [];
        data_ICS(i).raw(1:25) = [];
        data_ICS(i).sig(1:25) = [];
        offset = mean(data_ICS(i).sig(50:round(450 / data_ICS(i).freq / 1e6)));
        data_ICS(i).sigint = cumtrapz(data_ICS(i).sig - offset) * data_ICS(i).freq;
    end

    if ((isempty(data_NI) == 1) && (isempty(data_NJ) == 1) && (isempty(data_NK) == 1) && (isempty(data_ICS) == 1))
        magdata = [];
        return
    end

    %% make probe position array
    data_full = [data_NI; data_NJ; data_NK; data_ICS];

    if strcmp(date, '171018') || strcmp(date, '171019') || strcmp(date, '171020') || strcmp(date, '171023')
        data_full = [data_NI; data_NK; data_ICS];
    end

    if (strcmp(probe, 'NULL') == 1)
        Ntnum = 0;
        Nznum = 0;
        Ntarray = zeros(64, 1);
        Nzarray = zeros(64, 1);

        for i = 1:64

            for j = 1:length(data_full)
                s1 = ['Nt' num2str(i, '%02d')];
                s2 = ['Nz' num2str(i, '%02d')];

                if (strcmp(data_full(j).name, s1) == 1)
                    Ntarray(i) = j;
                    Ntnum = Ntnum + 1;
                elseif (strcmp(data_full(j).name, s2) == 1)
                    Nzarray(i) = j;
                    Nznum = Nznum + 1;
                end

            end

        end

        if ((Ntnum < 64) || (Nznum < 64))
            magdata = [];
            return
        end

    elseif (strcmp(probe, 'MID') == 1)
        Mtnum = 0;
        Mznum = 0;
        Mtarray = zeros(81, 1);
        Mzarray = zeros(81, 1);

        for i = 1:81

            for j = 1:length(data_full)
                s1 = ['Mt' num2str(i, '%02d')];
                s2 = ['Mz' num2str(i, '%02d')];

                if (strcmp(data_full(j).name, s1) == 1)
                    Mtarray(i) = j;
                    Mtnum = Mtnum + 1;
                elseif (strcmp(data_full(j).name, s2) == 1)
                    Mzarray(i) = j;
                    Mznum = Mznum + 1;
                end

            end

        end

        if ((Mtnum < 81) || (Mznum < 81))
            magdata = [];
            return
        end

    end

    if (exist([dir 'cdf\' num2str(date, '%06d') '\angle2.txt']) > 0)
        angle = load([dir 'cdf\' num2str(date, '%06d') '\angle2.txt']);
    elseif (exist([dir 'cdf\' num2str(date, '%06d') '\angle.txt']) > 0)
        angle = load([dir 'cdf\' num2str(date, '%06d') '\angle.txt']);
    else
        angle = zeros(81 + 64, 2);
        angle(:, 1) = linspace(1, 81 + 64, 81 + 64)';
        %disp('No angle.txt file!');
        handles.message.ForegroundColor = [1 0 0];
        handles.message.String = 'No angle.txt file   ';
    end

    calib = load('localdata\RBtcalib.mat');

    switch probe
        case 'NULL'
            %% NULL probe
            r = [0.11 0.17 0.23 0.29 0.35 0.41 0.47 ...
                    0.11 0.16 0.21 0.26 0.31 0.36 0.41 0.46 ...
                    0.11 0.17 0.23 0.29 0.35 0.41 0.47 ...
                    0.11 0.17 0.23 0.29 0.35 0.41 0.47 ...
                    0.11 0.17 0.23 0.29 0.35 0.41 0.47 ...
                    0.11 0.17 0.23 0.29 0.35 0.41 0.47 ...
                    0.11 0.17 0.23 0.29 0.35 0.41 0.47 ...
                    0.11 0.17 0.23 0.29 0.35 0.41 0.47 ...
                    0.11 0.17 0.23 0.29 0.35 0.41 0.47];
            z = [0.948 0.948 0.948 0.948 0.948 0.948 0.948 ...
                    0.873 0.873 0.873 0.873 0.873 0.873 0.873 0.873 ...
                    0.798 0.798 0.798 0.798 0.798 0.798 0.798 ...
                    0.723 0.723 0.723 0.723 0.723 0.723 0.723 ...
                    0.648 0.648 0.648 0.648 0.648 0.648 0.648 ...
                    0.573 0.573 0.573 0.573 0.573 0.573 0.573 ...
                    0.498 0.498 0.498 0.498 0.498 0.498 0.498 ...
                    0.423 0.423 0.423 0.423 0.423 0.423 0.423 ...
                    0.348 0.348 0.348 0.348 0.348 0.348 0.348];
            magdata_s = struct('Btm', cell(1), 'Btd', cell(1), 'Bttime', cell(1), 'Bzm', cell(1), 'Bzd', cell(1), 'Bztime', cell(1), 'angle', cell(1), ...
                'r', cell(1), 'z', cell(1), 'type', cell(1), 'Bz', cell(1), 'Bt', cell(1), 'Bzdev', cell(1), 'Btdev', cell(1), 'Bt2', cell(1));
            magdata = repmat(magdata_s, 64, 1);

            for i = 1:64
                magdata(i).Btm = data_full(Ntarray(i)).sigint;
                magdata(i).Btd = data_full(Ntarray(i)).sig;
                magdata(i).Bttime = data_full(Ntarray(i)).time;
                magdata(i).Bzm = data_full(Nzarray(i)).sigint;
                magdata(i).Bzd = data_full(Nzarray(i)).sig;
                magdata(i).Bztime = data_full(Nzarray(i)).time;
                magdata(i).angle = angle(i, 2) / 180 * pi;
                magdata(i).r = r(i);
                magdata(i).z = z(i);
                magdata(i).type = 'NULL';
            end

            for i = 1:64
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

        case 'MID'
            %% MID probe
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
            magdata_s = struct('Btm', cell(1), 'Btd', cell(1), 'Bttime', cell(1), 'Bzm', cell(1), 'Bzd', cell(1), 'Bztime', cell(1), 'anle', cell(1), ...
                'r', cell(1), 'z', cell(1), 'type', cell(1), 'Bz', cell(1), 'Bt', cell(1), 'Bzdev', cell(1), 'Btdev', cell(1), 'Bt2', cell(1));
            magdata = repmat(magdata_s, 81, 1);

            for i = 1:81
                magdata(i).Btm = data_full(Mtarray(i)).sigint;
                magdata(i).Btd = data_full(Mtarray(i)).sig;
                magdata(i).Bttime = data_full(Mtarray(i)).time;
                magdata(i).Bzm = data_full(Mzarray(i)).sigint;
                magdata(i).Bzd = data_full(Mzarray(i)).sig;
                magdata(i).Bztime = data_full(Mzarray(i)).time;

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

    end
