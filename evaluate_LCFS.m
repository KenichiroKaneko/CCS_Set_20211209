function err_LCFS = evaluate_LCFS(psi, REF, PARAM, CONFIG, CCR, CCZ, roop)
    % 磁気軸からLCFSまでの距離から誤差評価する関数
    % 放射状に線を引くことで距離を取得する
    % 何度もこの関数を呼び出す場合、その度にreferenceの計算をするのがリソースの無駄になるので
    % roopが2以上だったら保存したreferenceの値を読み込んで利用する

    % l 直線の長さ/s 直線を分割する数/ignore 内側から無視する数（中心付近は磁束が乱れる）/n 放射状に線を引く数
    l = 0.6;
    s = 1000;
    ignore = round(s * 0.1);
    n = 32;

    REF.Flux(:,1:5) = 1;
    REF01 = REF.Flux<0;


    if (roop > 1 && exist([PARAM.temporary_file_directory '/evaluate']))
        load([PARAM.temporary_file_directory '/evaluate'])
    else
        % ReferenceのLCFSの評価
        % 合体後プラズマ
        if PARAM.CCS == 1
            if CONFIG.ShowFig == 1
                figure()
                hold on
            end
            % 磁気軸
            [m r0] = min(min(REF.Flux));
            [m z0] = min(REF.Flux(:, r0));
            r0 = REF.R(r0);
            z0 = REF.Z(z0);
            for i = 1:n
                theta = 2 * pi / n * (i - 1) + pi / n;
                r = l * sin(theta); z = l * cos(theta);
                line_r = [r0, r0 + r]; line_z = [z0, z0 + z];
                % ある直線が作るpsiの断面を取得/
                [cx_ref, cy_ref, c_ref] = improfile(REF.R, REF.Z, REF01, line_r, line_z, s, 'bilinear');
                % LCFS(psi=0)との交点を探索
                [m, J] = min(abs(c_ref));
                I_ref(i) = J;
                if CONFIG.ShowFig == 1
                    % plot3(cx_ref, cy_ref, c_ref);
                    plot(cx_ref(I_ref(i)), cy_ref(I_ref(i)), 'm*')
                end
            end
            save([PARAM.temporary_file_directory '/evaluate'], 'I_ref')
        elseif PARAM.CCS == 2
            zhalf = round(length(REF.Z) / 2);
            if CONFIG.ShowFig == 1
                figure()
                hold on
            end
            for c = 1:2
                if c == 1
                    range = 1:zhalf;
                    zplus = 0;
                else
                    range = zhalf+1:length(REF.Z);
                    zplus = zhalf;
                end
                % 磁気軸
                [m r0] = min(min(REF.Flux(range, :)));
                [m z0] = min(REF.Flux(range, r0));
                R0(c) = REF.R(r0);
                Z0(c) = REF.Z(z0 + zplus);
                for i = 1:n
                    theta = 2 * pi / n * (i - 1) + pi / n;
                    r = l * sin(theta); z = l * cos(theta);
                    line_r = [R0(c), R0(c) + r]; line_z = [Z0(c), Z0(c) + z];
                    % ある直線が作るpsiの断面を取得 
                    [cx_ref, cy_ref, c_ref] = improfile(REF.R, REF.Z(range), REF01(range, :), line_r, line_z, s, 'bilinear');
                    % LCFS(psi=0)との交点を探索
                    [m, J] = min(abs(c_ref));
                    I_ref(c, i) = J;
                    if CONFIG.ShowFig == 1
                        plot3(cx_ref, cy_ref, c_ref);
                        plot(cx_ref(I_ref(c, i)), cy_ref(I_ref(c, i)), 'o')
                    end
                end
            end
            I_ref = reshape(I_ref, [1, 2*n]);
            save([PARAM.temporary_file_directory '/evaluate'], 'I_ref')
        end
    end
    
    psi(:,1:2) = 1;
    psi01 = psi<0;
    se = strel('disk',5);
    psi01 = imclose(psi01, se);

    if CONFIG.ShowFig == 1
        v = linspace(-20, 20, 101);
        % figure('Name', 'evaluateLCFS')
        % hold on
        % contour(CCR, CCZ, psi*1000, v)
        contour(CCR, CCZ, psi, [0 0], 'k', 'LineWidth', 2); % LCFS
        contour(REF.R, REF.Z, REF.Flux, [0 0], 'LineColor', 'm', 'LineWidth', 2);
    end

    if PARAM.CCS == 1
        for i = 1:n
            theta = 2 * pi / n * (i - 1) + pi / n;
            r = l * sin(theta); z = l * cos(theta);
            line_r = [r0, r0 + r]; line_z = [z0, z0 + z];
            % ある直線が作るpsiの断面を取得
            [cx,cy,c] = improfile(CCR, CCZ, psi01, line_r, line_z, s, 'bilinear');
            % LCFS(psi=0)との交点を探索
            [m, J] = min(abs(c));
            I(i) = J;
            if CONFIG.ShowFig
                % plot([line_r], [line_z], 'r');
                % plot(cx_ref(I_ref(i)), cy_ref(I_ref(i)), 'm*')
                plot(cx(I(i)), cy(I(i)), 'co');
            end
        end
    elseif PARAM.CCS == 2
        zhalf = round(size(psi, 1) / 2);
        for c = 1:2
            if c == 1
                range = 1:zhalf;
            else
                range = zhalf+1:size(psi, 1);
            end
        
            for i = 1:n
                theta = 2 * pi / n * (i - 1) + pi / n;
                r = l * sin(theta); z = l * cos(theta);
                line_r = [R0(c), R0(c) + r]; line_z = [Z0(c), Z0(c) + z];
                % ある直線が作るpsiの断面を取得
                [cx,cy,cc] = improfile(CCR, CCZ(range), psi01(range, :), line_r, line_z, s, 'bilinear');
                cc(1:ignore) = 0.01;
                [m, J] = min(abs(cc));
                I(c, i) = J;
                if CONFIG.ShowFig
                    plot([line_r], [line_z], 'r');
                    plot(cx(I(c, i)), cy(I(c, i)), 'o')
                end
            end
        end
        I= reshape(I, [1, 2*n]);
    end
    
    % matlabは1 origin なのでインデックスを-1する
    % 正しい距離を計算するのも、インデックスで計算するのも結果は同じ
    % LCFSの誤差評価
    err_LCFS = 1 / n * sum(abs(I_ref - I) ./ (I_ref-1)) * 100;

    figure()
    imshow(psi01)
end