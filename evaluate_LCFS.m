function err_LCFS = evaluate_LCFS(psi, REF, PARAM, CONFIG, CCR, CCZ, roop)
    % 磁気軸からLCFSまでの距離から誤差評価する関数
    % 放射状に線を引くことで距離を取得する
    % 何度もこの関数を呼び出す場合、その度にreferenceの計算をするのがリソースの無駄になるので
    % roopが2以上だったら保存したreferenceの値を読み込んで利用する

    % 磁気軸
    [m r0] = min(min(REF.Flux));
    [m z0] = min(REF.Flux(:, r0));
    r0 = REF.R(r0);
    z0 = REF.Z(z0);
    % l 直線の長さ/s 直線を分割する数/ignore 内側から無視する数（中心付近は磁束が乱れる）/n 放射状に線を引く数
    l = 0.6;
    s = 1000;
    ignore = round(s * 0.3);
    n = 32;

    if (roop > 1 && exist([PARAM.temporary_file_directory '/evaluate']))
        load([PARAM.temporary_file_directory '/evaluate'])
    else
        % ReferenceのLCFSの評価
        for i = 1:n
            theta = 2 * pi / n * (i - 1) + pi / n;
            r = l * sin(theta); z = l * cos(theta);
            line_r = [r0, r0 + r]; line_z = [z0, z0 + z];
            % ある直線が作るpsiの断面を取得
            [cx_ref, cy_ref, c_ref] = improfile(REF.R, REF.Z, REF.Flux, line_r, line_z, s, 'bilinear');
            c_ref(1:ignore) = [];
            % LCFS(psi=0)との交点を探索
            [m, J] = min(abs(c_ref));
            I_ref(i) = J + ignore;
        end
        save([PARAM.temporary_file_directory '/evaluate'], 'I_ref')
    end
    
    if CONFIG.ShowFig == 1
        v = linspace(-20, 20, 101);
        figure('Name', 'evaluateLCFS')
        hold on
        contour(CCR, CCZ, psi*1000, v)
        contour(CCR, CCZ, psi, [0 0], 'k', 'LineWidth', 2); % LCFS
        contour(REF.R, REF.Z, REF.Flux, [0 0], 'LineColor', 'm', 'LineWidth', 2);
    end
    for i = 1:n
        theta = 2 * pi / n * (i - 1) + pi / n;
        r = l * sin(theta); z = l * cos(theta);
        line_r = [r0, r0 + r]; line_z = [z0, z0 + z];
        % ある直線が作るpsiの断面を取得
        [cx,cy,c] = improfile(CCR, CCZ, psi, line_r, line_z, s, 'bilinear');
        c(1:ignore) = [];
        % LCFS(psi=0)との交点を探索
        [m, J] = min(abs(c));
        I(i) = J + ignore;
        if CONFIG.ShowFig
            plot([line_r], [line_z], 'r');
        end
    end
    
    % matlabは1 origin なのでインデックスを-1する
    % 正しい距離を計算するのも、インデックスで計算するのも結果は同じ
    % LCFSの誤差評価
    err_LCFS = 1 / n * sum(abs(I_ref - I) ./ (I_ref-1)) * 100;

end