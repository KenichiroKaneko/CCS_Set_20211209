%% loadwalldata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadwalldata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadwalldata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadwalldata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function WALL = loadwalldata(PARAM, CONFIG)

    %% Here, the geometory should be derived from input file
    %     R = [0.694D0, 0.694D0, 0.5985D0, 0.5985D0, 0.10815D0, 0.10815D0,...
    %         0.10815D0, 0.5985D0, 0.5985D0, 0.694D0,  0.694D0];
    %     Z = [0.0D0, 0.285D0, 0.285D0, 0.9985D0, 0.9985D0, 0.0D0,...
    %         -0.9985D0, -0.9985D0, -0.285D0, -0.285D0 , 0.0D0];

    %% CONDSHL_UTST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % 導体上の渦電流節点位置データの設定と読み取り
    %function [KNE,KNN,REV,ZEV,RES,ZES,KSE,KSN,REVN,ZEVN] = CONDSHL_UTST(WAHAHA,Nedp,NONC)

    MAXM = 13;
    Z = 0.15;
    INBOAD = 0.5;

    MSEC = [1, 1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1, 1];

    RSEC(1) = 0.694; ZSEC(1) = 0.0;
    RSEC(2) = 0.694; ZSEC(2) = 0.285;
    RSEC(3) = 0.5985; ZSEC(3) = 0.285;
    RSEC(4) = 0.5985; ZSEC(4) = 0.285 + Z;
    RSEC(5) = 0.5985; ZSEC(5) = 0.9985;
    RSEC(6) = 0.10815; ZSEC(6) = 0.9985;
    RSEC(7) = 0.10815; ZSEC(7) = INBOAD;
    RSEC(8) = 0.10815; ZSEC(8) = -INBOAD;
    RSEC(9) = 0.10815; ZSEC(9) = -0.9985;
    RSEC(10) = 0.5985; ZSEC(10) = -0.9985;
    RSEC(11) = 0.5985; ZSEC(11) = -0.285 - Z;
    RSEC(12) = 0.5985; ZSEC(12) = -0.285;
    RSEC(13) = 0.694; ZSEC(13) = -0.285;
    RSEC(14) = 0.694; ZSEC(14) = 0.0;
    %% Modified_MI 20211103
    % RSEC(14) = []; ZSEC(14) = [];
    % RSEC(1) = []; ZSEC(1) = [];
    MSEC = [1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1]; %original
    % MSEC = [1, 1, 2, 2, 2, 4, 2, 2, 2, 1, 1]; % ほどほどにノード点を増やした
    % MAXM = 11;
    %% Modified_MI 20211103 kokomade

    if 1 %CONFIG.DataType == 'sol'
        % 数値解は完全導体を仮定しているので外側の渦電流ノードを消さない
        % と思いきや凸部の渦電流を無視できそうなので全部else側に行かせてる
        % kaneko
        RSEC_list = [0.694, 0.694, 0.5985, 0.5985, 0.5985, 0.10185, 0.10185];
        RSEC_list = [RSEC_list fliplr(RSEC_list)];
        ZSEC_list = [0.0, 0.285, 0.285, 0.435, 0.9985, 0.9985, 0.6];
        ZSEC_list = [ZSEC_list -fliplr(ZSEC_list)];
        MAXM = length(RSEC_list) - 1;
        % MSEC = [1, 2, 1, 2, 3, 1, 2, 1, 3, 2, 1, 2, 1];
        MSEC = [2, 1, 1, 2, 2, 1, 3, 1, 2, 2, 1, 1, 2];
        MSEC = ones(MAXM);

        for i = 1:MAXM + 1
            RSEC(i) = RSEC_list(i);
            ZSEC(i) = ZSEC_list(i);
        end

        RWALL_list = [0.694, 0.694, 0.5985, 0.5985, 0.10185];
        ZWALL_list = [0, 0.285, 0.285, 0.9985, 0.9985];

    else
        % UTSTは外側で渦電流が流れにくいので渦電流ノードを一部消す
        % kaneko
        % 凸壁なし
        RSEC_list = [0.694, 0.5985, 0.5985, 0.5985, 0.10185, 0.10185];
        ZSEC_list = [0.285, 0.285, 0.435, 0.9985, 0.9985, 0.6];
        % MSEC = [1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1];
        MSEC = [1, 1, 2, 3, 1, 3, 1, 3, 2, 1, 1];
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

        RWALL_list = [0.694, 0.5985, 0.5985, 0.10185];
        ZWALL_list = [0.285, 0.285, 0.9985, 0.9985];

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

    WALL.KNE = 0;

    for I = 1:MAXM
        WALL.KNE = WALL.KNE + MSEC(I);
        CM = 2 * MSEC(I);
        DELR = (RSEC(I + 1) - RSEC(I)) / CM;
        DELZ = (ZSEC(I + 1) - ZSEC(I)) / CM;

        for J = 1:2 * MSEC(I)
            II = II + 1;
            WALL.REV(II) = WALL.REV(II - 1) + DELR;
            WALL.ZEV(II) = WALL.ZEV(II - 1) + DELZ;
        end

    end

    KNM = II - 1;

    RWALL_list = [RWALL_list(end) RWALL_list];
    ZWALL_list = [ZWALL_list(end) ZWALL_list];

    %% ！！！！！！！　 非適合要素節点計算の予定地　　　！！！！！！！
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
    %***
    %cd      KSE=1
    %prompt = '電流シート上に境界要素を配置する？ (Yes/No)=(1/0)\n';
    %KSE = input(prompt);

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
