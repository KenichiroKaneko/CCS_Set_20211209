function WALL = loadwalldata2(PARAM, CONFIG)


    if CONFIG.DataType == 'sol'
        % 数値解は完全導体を仮定しているので外側の渦電流ノードを消さない
        % kaneko
        RSEC_list = [0.694, 0.694, 0.5985, 0.5985, 0.5985, 0.10185, 0.10185];
        RSEC_list = [RSEC_list fliplr(RSEC_list)];
        ZSEC_list = [0.0, 0.285, 0.285, 0.435, 0.9985, 0.9985, 0.6];
        ZSEC_list = [ZSEC_list -fliplr(ZSEC_list)];
        MAXM = length(RSEC_list) - 1;
        % MSEC = [1, 2, 1, 2, 3, 1, 2, 1, 3, 2, 1, 2, 1];
        % MSEC = [2, 1, 1, 2, 2, 1, 3, 1, 2, 2, 1, 1, 2];
        MSEC = [1, 1, 1, 2, 2, 1, 3, 1, 2, 2, 1, 1, 1];
        % MSEC = [1, 2, 1, 1, 1, 1, 3, 1, 1, 1, 1, 2, 1];
        % MSEC = [0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 2, 0];
        % MSEC = ones(MAXM);



        % 配置1
        % RSEC_list = [0.69, 0.69, 0.59, 0.59, 0.10];
        % RSEC_list = [RSEC_list, 0.10, fliplr(RSEC_list)];
        % ZSEC_list = [0.00, 0.285, 0.285, 1.0, 1.0];
        % ZSEC_list = [ZSEC_list, 0, -fliplr(ZSEC_list)];
        % MAXM = length(RSEC_list) - 1;
        % MSEC = [PARAM.eddyNodes fliplr(PARAM.eddyNodes)];

        % 配置2
        RSEC_list = [0.694, 0.694, 0.5985, 0.5985, 0.5985, 0.10185, 0.10185];
        RSEC_list = [RSEC_list, fliplr(RSEC_list)];
        ZSEC_list = [0.0, 0.285, 0.285, 0.435, 0.9985, 0.9985, 0.6];
        ZSEC_list = [ZSEC_list, -fliplr(ZSEC_list)];
        MAXM = length(RSEC_list) - 1;
        MSEC = [PARAM.eddyNodes fliplr(PARAM.eddyNodes(1:end-1))];

        for i = 1:MAXM + 1
            RSEC(i) = RSEC_list(i);
            ZSEC(i) = ZSEC_list(i);
        end

        RWALL_list = [0.694, 0.694, 0.5985, 0.5985, 0.10185];
        ZWALL_list = [0, 0.285, 0.285, 0.9985, 0.9985];
        % RWALL_list = [0.69, 0.59, 0.59, 0.10];
        % ZWALL_list = [0.3, 0.3, 1, 1];

    else
        % UTSTは外側で渦電流が流れにくいので渦電流ノードを一部消す
        % kaneko
        % 凸壁後で消去
        RSEC_list = [0.694, 0.694, 0.5985, 0.5985, 0.5985, 0.10185, 0.10185];
        ZSEC_list = [0.0, 0.285, 0.285, 0.435, 0.9985, 0.9985, 0.6];
        MSEC = [0, 2, 2, 2, 2, 1, 3, 1, 2, 2, 2, 2, 0];
        % 凸壁なし
        % RSEC_list = [0.694, 0.5985, 0.5985, 0.5985, 0.10185, 0.10185];
        % ZSEC_list = [0.285, 0.285, 0.435, 0.9985, 0.9985, 0.6];
        % MSEC = [1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1];
        % MSEC = [1, 1, 2, 3, 1, 3, 1, 3, 2, 1, 1];
        % 1 PF部分削減
        % RSEC_list = [0.694, 0.5985, 0.5985, 0.10185, 0.10185];
        % ZSEC_list = [0.285, 0.285, 0.9985, 0.9985, 0.6];
        % MSEC = [1, 3, 2, 1, 3, 1, 2, 3, 1];
        % 2 凸なし
        % RSEC_list = [0.5985, 0.5985, 0.10185, 0.10185];
        % ZSEC_list = [0.285, 0.9985, 0.9985, 0.6];
        % MSEC = [4, 2, 1, 3, 1, 2, 4];
        % 3 内側削減
        % RSEC_list = [0.694, 0.5985, 0.5985, 0.10185];
        % ZSEC_list = [0.285, 0.285, 0.9985, 0.9985];
        % MSEC = [1, 1, 1, 2, 1, 1, 1];
        % MSEC = [2, 3, 3, 4, 3, 3, 2];
        % 4 凸なし内側削減
        % RSEC_list = [0.5985, 0.5985, 0.10185];
        % ZSEC_list = [0.285, 0.9985, 0.9985];
        % MSEC = [2, 1, 3, 1, 2];
        % MSEC = [2, 3, 3, 4, 3, 3, 2];

        RSEC_list = [RSEC_list fliplr(RSEC_list)];
        ZSEC_list = [ZSEC_list -fliplr(ZSEC_list)];
        % RSEC_list = [RSEC_list 0.10185 fliplr(RSEC_list)];
        % ZSEC_list = [ZSEC_list 0 -fliplr(ZSEC_list)];
        MAXM = length(RSEC_list) - 1;
        % MSEC = ones(MAXM) * 3;

        for i = 1:MAXM + 1
            RSEC(i) = RSEC_list(i);
            ZSEC(i) = ZSEC_list(i);
        end

        % RWALL_list = [0.694, 0.5985, 0.5985, 0.10185];
        % ZWALL_list = [0.285, 0.285, 0.9985, 0.9985];
        RWALL_list = [0.69, 0.59, 0.59, 0.10];
        ZWALL_list = [0.3, 0.3, 1, 1];

    end

    %%
    % eddycurrentのplot用にノード点の距離を取得
    RWALL_list = [RWALL_list fliplr(RWALL_list)];
    ZWALL_list = [ZWALL_list -fliplr(ZWALL_list)];
    RWALL_list2 = abs(RWALL_list(1:end - 1) - RWALL_list(2:end));
    ZWALL_list2 = abs(ZWALL_list(1:end - 1) - ZWALL_list(2:end));
    WALL_dist_list = [];
    WALL.RSEC = RSEC;
    WALL.ZSEC = ZSEC;


    for i = 1:length(RWALL_list) - 1
        WALL_dist_list(i) = sum(RWALL_list2(1:i)) + sum(ZWALL_list2(1:i));
    end

    WALL.WALL_dist_list = [0 WALL_dist_list];

    RWALL_list = [RWALL_list(end) RWALL_list];
    ZWALL_list = [ZWALL_list(end) ZWALL_list];
    WALL.RWALL_list = RWALL_list;
    WALL.ZWALL_list = ZWALL_list;
    % plot用のコードここまで

    %%

    % 真空容器上の渦電流節点の設定
    % KNE=真空容器の分割境界要素数,  KNM=真空容器上のメッシュ点数,  KNN=真空容器上の節点数

    II = 1;
    WALL.REV(II) = RSEC(1);
    WALL.ZEV(II) = ZSEC(1);
    WALL.KNE = sum(MSEC);

    for I = 1:MAXM
        if MSEC(I) == 0
            DELR = (RSEC(I + 1) - RSEC(I));
            DELZ = (ZSEC(I + 1) - ZSEC(I));
            II = II + 1;
            WALL.REV(II) = WALL.REV(II - 1) + DELR;
            WALL.ZEV(II) = WALL.ZEV(II - 1) + DELZ;
        else
            CM = 2 * MSEC(I);
            DELR = (RSEC(I + 1) - RSEC(I)) / CM;
            DELZ = (ZSEC(I + 1) - ZSEC(I)) / CM;
            for J = 1:2 * MSEC(I)
                II = II + 1;
                WALL.REV(II) = WALL.REV(II - 1) + DELR;
                WALL.ZEV(II) = WALL.ZEV(II - 1) + DELZ;
            end
        end
    end

    if (MSEC(1) + MSEC(end))==0
        % MSECの最初と最後が両方０なら点を削除
        WALL.REV(end) = [];
        WALL.ZEV(end) = [];
        WALL.KNE = WALL.KNE - 1;
        WALL.REV(1) = [];
        WALL.ZEV(1) = [];
        WALL.KNE = WALL.KNE - 1;
    else
        % 最初と最後の点を同じにして繋げる
        WALL.REV(end) = WALL.REV(1);
        WALL.ZEV(end) = WALL.ZEV(1);
    end
    
    KNM = length(WALL.REV);

    RWALL_list = [RWALL_list(end) RWALL_list];
    ZWALL_list = [ZWALL_list(end) ZWALL_list];

    %% ！非適合要素節点計算の予定地！
    % 真空容器上の非適合要素節点座標の作成
    if (PARAM.NONC > 0)
        WALL.KNN = WALL.KNE * 3;
        I = 1:WALL.KNE;
        WALL.REVN(3 * I - 2) = (5.D0 * WALL.REV(2 * I - 1) + 5.D0 .* WALL.REV(2 * I) - WALL.REV(2 * I + 1)) / 9.D0;
        WALL.ZEVN(3 * I - 2) = (5.D0 * WALL.ZEV(2 * I - 1) + 5.D0 .* WALL.ZEV(2 * I) - WALL.ZEV(2 * I + 1)) / 9.D0;
        WALL.REVN(3 * I - 1) = WALL.REV(2 * I);
        WALL.ZEVN(3 * I - 1) = WALL.ZEV(2 * I);
        WALL.REVN(3 * I) = (5.D0 * WALL.REV(2 * I + 1) + 5.D0 .* WALL.REV(2 * I) - WALL.REV(2 * I - 1)) / 9.D0;
        WALL.ZEVN(3 * I) = (5.D0 * WALL.ZEV(2 * I + 1) + 5.D0 .* WALL.ZEV(2 * I) - WALL.ZEV(2 * I - 1)) / 9.D0;
    else
        WALL.KNN = KNM;
    end

    % 安定化板上の渦電流節点の設定
    % KSE=安定化板の分割境界要素数,  KSN=安定化板上の節点数

    WALL.KSE = 0;
    % centerR = sum(PARAM.R0) / length(PARAM.R0);
    % centerZ = sum(PARAM.Z0) / length(PARAM.Z0);
    % WALL.RES = [centerR - 0.05, centerR, centerR + 0.05];
    % WALL.ZES = [centerZ, centerZ, centerZ];

    %***
    WALL.KSN = WALL.KSE * 2 + 1;

    if (WALL.KSE == 0)
        WALL.KSN = 0;
    end

    if (WALL.KSE > 0)
        fp = fopen([PARAM.temporary_file_directory '/StabilizerPoints.txt'], 'w'); % 114
        WALL.RES = [centerR - 0.1, centerR - 0.05, centerR, centerR + 0.05, centerR + 0.1];
        WALL.ZES = centerZ * ones(length(WALL.RES));
        % WALL.ZES = zeros(length(WALL.RES));
        % WALL.RES(1) = 0.25;
        % WALL.RES(2) = 0.28;
        % WALL.RES(3) = 0.31;
        % WALL.ZES(1) = 0;
        % WALL.ZES(2) = 0;
        % WALL.ZES(3) = 0;

        for I = 1:WALL.KSN
            fprintf(fp, '%d %d\n', WALL.RES(I), WALL.ZES(I));
        end

        fclose(fp);
    else
        WALL.RES = [];
        WALL.ZES = [];
    end

end

%% loadwalldata kokomade!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
