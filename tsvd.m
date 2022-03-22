function FFout = tsvd(PARAM, CONFIG, A, FF, FC, ...
    SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL)
    % 打ち切り特異値分解

    LTikh = 0;
    GAM0 = 0;

    [M, N] = size(A); % 未知数、既知数
    if (M < N)
        disp([M N])
        error('You must augment A with extra zero rows.', N)
    end
    if CONFIG.ShowFig
        disp(['既知数 : ' num2str(M) '   未知数 : ' num2str(N)])
    end
    NAPB = SENSOR_TPRB.NUM + SENSOR_NPRB.NUM;
    NFLX = SENSOR_FLXLP.NUM;
    NCCN = CCSDAT.NCCN;
    KNN = WALL.KNN;
    KSN = WALL.KSN;

    % A = U*S*V'
    [U, S, V] = svd(A);
    W = diag(S);
    KUP_L = Lcurve(PARAM, CONFIG, A, W, V, U, FF);

    % PARAMでKUPを指定していなかったらL-curve法で指定する
    if PARAM.KUP == 0
        PARAM.KUP = KUP_L;
        KUP = KUP_L;
    else
        KUP = PARAM.KUP;
    end

    % KUP = max(KUP, 75);

    if CONFIG.ShowFig
        disp(['Lcurve : ', num2str(KUP_L) ', KUP : ', num2str(KUP)])
    end
    % 特異値行列の逆行列を計算 S(MxN)→S*(NxM)
    W_inv = zeros(N, M);

    for i = 1:KUP;
        W_inv(i, i) = 1 / W(i);
    end

    % 求めるベクトル
    P_solved = V * W_inv * inv(U) * FF';
    X = P_solved;
    FFout = X;

    if CONFIG.ShowFig2
        cal_MRE(PARAM, CONFIG, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, A, S, V, U, X, FC, FF)
    end
end