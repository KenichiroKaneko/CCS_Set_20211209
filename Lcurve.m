% Lcurve法で打切り特異値分解の打切り項数を決定する

function KUP = Lcurve(PARAM, CONFIG, A, W, V, U, B)
    % A：特異値分解される行列
    % W：特異値の行列
    % V：[U S V] = svd(A)のV
    % U：[U S V] = svd(A)のU
    % X：求める未知変数のベクトル
    % B：Ap = q の q
    % M：既知の数、N：未知の数
    [M N] = size(A);

    B = B';
    % W = diag(W); % 対角成分のみの行列
    W = sort(W, 'descend');
    W_length = length(W);

    for i = 1:W_length

        % 特異値行列の逆行列を計算 S(MxN)→S*(NxM)
        W_inv = zeros(N, M);

        for j = 1:i
            W_inv(j, j) = 1 / W(j);
        end

        % 求めるベクトル
        P_solved = V * W_inv * inv(U) * B;

        % 残差のノルムの計算
        residual(i) = norm(A * P_solved - B).^2;

        % 求めたベクトルのノルム
        P_norm(i) = norm(P_solved).^2;
    end


    % figure('Name', '特異値')
    % semilogy(1:75, W)

    KUP = cal_curvature(residual, P_norm, CONFIG);
end

function KUP = cal_curvature(X, Y, CONFIG);

    x = log10(X);
    y = log10(Y);

    % 前後の無視する数/点群の数、num < cut
    low_cut = 40;
    high_cut = 0;
    num = 1;

    len = length(X);
    m = 0;
    KUP = 1;

    % 原点をずらす
    x1 = x(len);
    y1 = y(1);
    x = x - x1;
    y = y - y1;

    %% 曲率最大求める方法その１
    % y=-x-10と一番近い点 ax + by + c = 0
    y_intercept = 1;
    a = 1;
    b = 1;
    a = -(y(low_cut) - y(len-high_cut));
    b = x(low_cut) - x(len-high_cut);
    c = y_intercept;
    j = 1;
    lx = [0 -c/a * 10];
    ly = [-c/b * 10 0];

    for i = low_cut+1:length(x) - high_cut
        d(j) = abs(a * x(i) + b * y(i) + c) / norm([a b]);
        j = j + 1;
    end

    %% 曲率最大求める方法その２
    % % 原点と一番近い点
    % for i = low_cut+1:length(x) - high_cut
    %     d(j) = x(i)^2 + y(i)^2;
    % end
    % figure()
    % plot(d)

    [My Iy] = min(d);
    Iy = Iy + low_cut;
    KUPIndex_Y = len + 1 - Iy;

    % ｘーｙと曲率最大のポイント赤色
    % figure()
    % plot(x, y, 's')
    % hold on
    % scatter(x(Iy), y(Iy), "y*")
    % hold off

    KUP = Iy;


    if CONFIG.ShowFig

        figure()
        loglog(X, Y, 's')
        hold on
        % plot(x, y, 's')
        plot(X(KUP), Y(KUP), 'mo', 'MarkerSize', 8);
        title("L-curve");
        xlabel("||Dp^* - g||", 'FontWeight','bold')
        ylabel("||p^*||", 'FontWeight','bold')
        % axis equal

        % for i = 5:5:len
        %     % text(residual(i), (P_norm(i)), ['\leftarrow', num2str(i)])
        %     text(X(i), (Y(i)), ['\leftarrow', num2str(i)])
        % end

    end
    % error('error description', A1)
end
