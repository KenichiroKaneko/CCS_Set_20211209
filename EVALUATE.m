% 最小２乗誤差を返すプログラム
function MSE = EVALUATE(psi, REF, PARAM, CONFIG, CCSDAT, CCR, CCZ)
    % psi 再構成結果の磁束分布
    % REF referenceの情報（磁束分布、R座標、Z座標）
    % PARAM, CONFIG 設定ファイル
    % CCSDAT CCS面に関するデータ（中心位置や大きさなど）
    % CCR,CCZ 再構成結果のR座標、Z座標のリスト

    % limiterRef = 0.1656627;
    limiterRef = -10;
    limiterRecon = -10;
    CONFIG.ShowFig=1;
    PARAM.CCS = 1;
    if CONFIG.ShowFig
        figure()
        contour(CCR, CCZ, psi, '-k', 'LevelStep', 0.0003);
        xlim([0.1 0.9]);
        ylim([-1 1]);
        xlabel({'r (m)'});
        ylabel({'z (m)'});
        title("Reconstructed LCFS")
        hold on
        axis equal
    end


    if 0 %PARAM.CCS > 1 && abs(PARAM.Z0(1) - PARAM.Z0(2)) > 0.8
        % 2つのプラズマが離れている場合

        for i = 1:2
            j = i - 1;
            range1 = 1 + j * 50:51 + j * 50;
            range2 = 1 + j * 1016:1017 + j * 1016;

            % 再構成した磁束のLCFSを探索
            [LCFSrecon, disuseNum1, disuseNum2, LCFS_min_recon] = detectLCFS(PARAM, CONFIG, CCSDAT, psi(range1, :), CCR, CCZ(range1), limiterRecon, "reconstruction");
            % referenceのLCFSを探索
            [LCFSref, MAxisR, MAxisZ, LCFS_min_ref] = detectLCFS(PARAM, CONFIG, CCSDAT, REF.Flux(range2, :), REF.R, REF.Z(range2), limiterRef, "reference");
            % LCFSを描画
            % contour(CCR, CCZ, psi(range1, :), [LCFS_min_recon LCFS_min_recon], 'LineColor', 'c', 'LineWidth', 2);
            % contour(REF.R, REF.Z, REF.Flux(range2, :), [LCFS_min_ref LCFS_min_ref], 'LineColor', 'm', 'LineWidth', 2);
            % 求めた2つの磁気面上の何点かと磁気軸からの距離をそれぞれ計算し、２乗誤差を計算
            MSE(i) = calcError(PARAM, CONFIG, LCFSrecon, LCFSref, MAxisR, MAxisZ, LCFS_min_recon, LCFS_min_ref, CCR, CCZ(range1), REF.R, REF.Z(range2));
        end

    else
        % 2つのプラズマが接近した場合
        [LCFSrecon, disuseNum1, disuseNum2, LCFS_min_recon] = detectLCFS(PARAM, CONFIG, CCSDAT, psi, CCR, CCZ, limiterRecon, "reconstruction");
        [LCFSref, MAxisR, MAxisZ, LCFS_min_ref] = detectLCFS(PARAM, CONFIG, CCSDAT, REF.Flux, REF.R, REF.Z, limiterRef, "reference");
        % LCFSを描画
        % contour(CCR, CCZ, psi, [LCFS_min_recon LCFS_min_recon], 'LineColor', 'c', 'LineWidth', 2);
        % contour(REF.R, REF.Z, REF.Flux, [LCFS_min_ref LCFS_min_ref], 'LineColor', 'm', 'LineWidth', 2);
        % 求めた2つの磁気面上の何点かと磁気軸からの距離をそれぞれ計算し、２乗誤差を計算
        MSE = calcError(PARAM, CONFIG, LCFSrecon, LCFSref, MAxisR, MAxisZ, LCFS_min_recon, LCFS_min_ref, CCR, CCZ, REF.R, REF.Z);
    end

    MSE = sum(MSE) / length(MSE);

    if CONFIG.ShowFig
        legend("Recon flux", "Ref LCFS", "Recon LCFS");
        str = sprintf('MSE= %.3e', MSE);
        % text(0.4, -0.8, str, 'FontSize', 10, 'BackgroundColor', "white")
    end
    hold off
end

function [STRUCTURE, MAxisR, MAxisZ, LCFS_min] = detectLCFS(PARAM, CONFIG, CCSDAT, psi, psiR, psiZ, limiter, type)
    % PARAM,CONFIG 設定ファイル
    % CCSDAT CCS面に関するデータ（中心位置や大きさなど）
    % LCFSを探す対象の磁束分布
    % psiR,psiZ LCFSを探す対象の磁束分布のR座標、Z座標
    % limiter リミターに関する調整
    % type reference/reconstruction

    psi_temp = psi;
    psi_inside = psi;
    plottype = 0;

    % 磁気軸の中心を探す→refでは利用されるがreconでは要らない
    [M1 MAxisR] = min(min(psi));
    [M2 MAxisZ] = min(psi(:, MAxisR));
    % scatter(psiR(MAxisR), psiZ(MAxisZ));

    if limiter == -1
        % リミターがない場合→０未満の磁束の最大値
        LCFS_min = max(max(psi_temp(psi_temp < 0))) - 0.0005;
    elseif limiter == -10
        % 再構成でもリミター位置での磁束を元に計算することにしたので後から追加した
        limiter = 0.1656627;
        [m limiterR] = min(abs(psiR - limiter));
        [m limiterZ] = min(psi(:, limiterR));
        LCFS_min = m;

    else
        % リミターがある場合→リミター位置の磁束
        % scatter(psiR(MAxisR), psiZ(MAxisZ), "*m");
        [m limiterR] = min(abs(psiR - limiter));
        LCFS_min = psi(MAxisZ, limiterR);
        % scatter(psiR(limiterR), psiZ(MAxisZ), "*")
    end

    if type == "reference"
        line = "-k";
    else
        line = "-r";
    end

    if plottype

        if type == "reconstruction" & CONFIG.ShowFig
            % contour(psiR, psiZ, psi, [LCFS_min LCFS_min], 'LineColor', 'r', 'LineWidth', 2)
        end

        STRUCTURE.nouse = 0;

    else

        psi_inside(psi_temp > LCFS_min) = 0; % LCFSの外側を０にする
        [Z, R] = ind2sub(size(psi_temp), find(psi_inside)); % LCFSの内側の座標

        if 0 %type == "reconstruction"
            % % CCNを広げるiranaikamo
            l = 1;
            n = length(CCSDAT.RCCN(1, :));

            for i = 1:n
                theta = 2 * pi / n * (i - 1) + pi / n;
                r = l;
                z = l;
                CCSDAT.RCCN(:, i) = CCSDAT.RCCN(:, i) + r;
                CCSDAT.ZCCN(:, i) = CCSDAT.ZCCN(:, i) + z;
            end

            if PARAM.CCS == 2
                plgn1 = polyshape(CCSDAT.RCCN(1, :), CCSDAT.ZCCN(1, :));
                plgn2 = polyshape(CCSDAT.RCCN(2, :), CCSDAT.ZCCN(2, :));
                plot(plgn1);
                plot(plgn2);
                scatter(psiR(R), psiZ(Z));
                in1 = inpolygon(psiR(R), psiZ(Z), CCSDAT.RCCN(1, :), CCSDAT.ZCCN(1, :));
                in2 = inpolygon(psiR(R), psiZ(Z), CCSDAT.RCCN(2, :), CCSDAT.ZCCN(2, :));
                in = in1 | in2;
            else
                % plgn1 = polyshape(CCSDAT.RCCN(1, :), CCSDAT.ZCCN(1, :));
                % plot(plgn1);
                % scatter(psiR(R), psiZ(Z));
                in = inpolygon(psiR(R), psiZ(Z), CCSDAT.RCCN(1, :), CCSDAT.ZCCN(1, :));
            end

            R = R(in);
            Z = Z(in);
        end

        j = boundary(R, Z, 0.1);
        STRUCTURE.j = j;
        STRUCTURE.R = R;
        STRUCTURE.Z = Z;
        % scatter(psiR(R), psiZ(Z), 'k');
    end

    % 他の境界プロット法
    % [B, L] = bwboundaries(psi_inside, 'noholes');

    % for k = 1:length(B)
    %     border = B{k};
    %     plot(psiR(border(:, 2)), psiZ(border(:, 1)), ':b', 'LineWidth', 1)
    % end
end

function MSE = calcError(PARAM, CONFIG, LCFSrecon, LCFSref, MAxisR, MAxisZ, LCFS_min_recon, LCFS_min_ref, CCR, CCZ, REFR, REFZ)
    % PARAM,CONFIG 設定ファイル
    % CCSDAT CCS面に関するデータ（中心位置や大きさなど）
    % LCFSrecon, LCFSref 構造体 jは２次元配列、境界の座標 R,Zは座標
    % MAxisR,MAxisZ referenceの磁気軸の中心 ここを基準に誤差を計測する
    % LCFS_min_recon, LCFS_min_ref LCFSを描く時の基準の値 この値の投稿線がLCFS
    % CCR,CCZ 再構成結果のR座標、Z座標のリスト
    % REFR, REFZ referencdのR座標、Z座標のリスト

    poly2 = polyshape(REFR(LCFSref.R(LCFSref.j)), REFZ(LCFSref.Z(LCFSref.j)));
    poly1 = polyshape(CCR(LCFSrecon.R(LCFSrecon.j)), CCZ(LCFSrecon.Z(LCFSrecon.j)));

    if CONFIG.ShowFig
        figure()
        hold on
        plot(poly2, 'FaceAlpha', 0.4, 'FaceColor', "m", 'EdgeColor', "m", "LineWidth", 3)
        plot(poly1, 'FaceAlpha', 0.4, 'FaceColor', "c", 'EdgeColor', "c", "LineWidth", 3)
        scatter(REFR(MAxisR), REFZ(MAxisZ));
        hold off
    end

    p = 32;
    % delta = 1.0e-8;
    delta = 0;
    l = 1;
    Z1 = REFZ(MAxisZ);
    R1 = REFR(MAxisR);

    for i = 1:p
        theta = 2 * (i + delta) * pi / p;

        % 磁気軸から結ぶ点
        R2 = R1 + l;
        Z2 = Z1 + l;
        % ２点を結ぶ直線
        lineseg = [R1 Z1; R2 Z2];

        % 再構成したLCFSと直線の交わり
        [in, out] = intersect(poly1, lineseg);
        x1 = in(1, 1);
        y1 = in(1, 2);
        x2 = in(end, 1);
        y2 = in(end, 2);
        L_recon = sqrt((x1 - x2)^2 + (y1 - y2)^2);
        % plot(in(:, 1), in(:, 2), 'b')
        % plot(out(:, 1), out(:, 2), 'r')
        % scatter(in(end, 1), in(end, 2), "ro")
        % scatter(in(1, 1), in(1, 2), "ro")

        % referenceLCFSと直線の交わり
        [in, out] = intersect(poly2, lineseg);
        x1 = in(1, 1); y1 = in(1, 2); x2 = in(end, 1); y2 = in(end, 2);
        L_ref = sqrt((x1 - x2)^2 + (y1 - y2)^2);
        % plot(in(:, 1), in(:, 2), 'b')
        % plot(out(:, 1), out(:, 2), 'r')
        % scatter(in(end, 1), in(end, 2), "yo")
        % scatter(in(1, 1), in(1, 2), "yo")

        error(i) = abs(L_recon - L_ref);
        L_ref_sum(i) = L_ref;
        % plot(x, y, "k")

    end

    MSE = sum(error) / p;
    MSE = sum(error) / sum(L_ref_sum);

end
