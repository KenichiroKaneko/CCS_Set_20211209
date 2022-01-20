%function [C,X,GET] = KUPCHK_matlab(PARAM,A,B,U,V,W,M,N,MP,NP,NAPB,NFLX,NCCN)
function [C, X, GET] = KUPCHK_matlab(PARAM, A, B, U, V, W, NAPB, NFLX, NCCN)
    [M, N] = size(A);
    C = zeros(1, N); %MP
    WTMP = zeros(1, 200);
    AIC0 = zeros(1, 200);
    AIC1 = zeros(1, 200);
    GET = zeros(1, N); %NP
    %
    %  ****  AX=B  ****
    %  Input==>  A(M,N): System Matrix,
    %            B(M): Inhomogeneous Term Vector
    %  Output==> X(N): Unknown Vector to be Solved
    %   M=Number of data points or equations,   N=Number of unknowns
    %   MP=Max. capacity of M,                  NP=Max. capacity of N
    % *****************************************************************
    %
    %
    % fid\d{3} = fopen([PARAM.temporary_file_directory '/@KUPCHK00.txt'], 'w'); %200
    % fid\d{3} = fopen([PARAM.temporary_file_directory '/@AIC0.txt'], 'w'); %201
    % fid\d{3} = fopen([PARAM.temporary_file_directory '/@AIC1.txt'], 'w'); %202
    % fid\d{3} = fopen([PARAM.temporary_file_directory '/@GETAvar.txt'], 'w'); %203
    %
    % fid\d{3} = fopen([PARAM.temporary_file_directory '/@NB_BNSN.txt'], 'w'); %204
    % fid\d{3} = fopen([PARAM.temporary_file_directory '/@NP_BNSN.txt'], 'w'); %205
    % fid\d{3} = fopen([PARAM.temporary_file_directory '/@NC_BNSN.txt'], 'w'); %206
    % fid\d{3} = fopen([PARAM.temporary_file_directory '/@NPARADD.txt'], 'w'); %207
    %
    AI0MIN = 1.0D10;
    AI1MIN = 1.0D10;

    for KUP = 1:N
        WTMP(1:N) = W(1:N);

        if (KUP ~= N)
            WTMP(KUP + 1:N) = 0.0D0;
        end

        %********** GETA iteration (START) *****************
        ITMX = 1 + PARAM.GETA_YN * 999;
        %ITMX = 1000;
        DELGE = 0.0D0;
        OLDEL = 1.0D0;
        GETA = 0.0D0;
        OLDGT = GETA;
        OMGA = 1.9D0;
        %           OMGA = 1.D0;
        for IT = 1:ITMX
            GETA = GETA + DELGE;
            GETA = OMGA * GETA + (1.0D0 - OMGA) * OLDGT;
            C(1:M) = B(1:M);
            C(NAPB + 1:NAPB + NFLX) = C(NAPB + 1:NAPB + NFLX) - GETA;
            %
            %               [X] = SVBKSB(U,WTMP,V,M,N,MP,NP,C); %OK
            [X] = SVBKSB(U, WTMP, V, C); %OK
            %
            %               BNSNALL = 0.0D0;
            %               BNSN_B = 0.0D0;
            %               BNSN_P = 0.0D0;
            %               BNSN_C = 0.0D0;
            %               DELGE = 0.0D0;
            %
            K = 1:M;
            J = 1:N;
            C(K) = A(K, J) * X(J)';
            E(K) = (C(K) - B(K) + GETA) .* and(K > NAPB, K <= (NFLX + NAPB)) + (C(K) - B(K)) .* or(K <= NAPB, K > (NFLX + NAPB)); % ! GETAを引くのはFlux loopのみ !!
            DELGE = sum((B(K) - C(K) - GETA) .* and(K > NAPB, K <= (NFLX + NAPB))); %  ! DELGE評価はFlux loopのみ
            BNSNALL = sum(E(K) .* E(K));
            BNSN_B = sum(E(K) .* E(K) .* (K <= NAPB));
            BNSN_P = sum(E(K) .* E(K) .* and(K > NAPB, K <= (NFLX + NAPB)));
            BNSN_C = sum(E(K) .* E(K) .* (K > (NFLX + NAPB)));
            CNT0 = NFLX;
            DELGE = DELGE / CNT0;
            EPSDEL = abs((DELGE - OLDEL) / OLDEL);

            if and(IT > 2, EPSDEL < 0.1e0)
                OMGA = OMGA + 5.0D0;
                % % fprintf(% fid\d{3}, '%s %d\r\n', '    ******  OMGA changed to', OMGA);
            end

            OLDEL = DELGE;
            OLDGT = GETA;
            EPS = abs(DELGE / GETA);
            % fprintf(% fid\d{3}, '%d %d %d %d\r\n', KUP, IT, GETA, DELGE);
            ITEND = IT;

            if (EPS < 1.0E-5)
                break
            end

        end

        % fprintf(% fid\d{3}, '%d %d %d %d\r\n', KUP, ITEND, GETA, DELGE);
        % fprintf(% fid\d{3}, '%d %d\r\n', KUP, GETA);
        GET(KUP) = GETA;
        %
        %********** GETA iteration (END) *******************
        %
        CN = N;
        CM = M;
        CNAPB = NAPB;
        CNFLX = NFLX;
        CNCCN = NCCN;
        BNSNALL = BNSNALL / CM;
        BNSN_B = BNSN_B / CNAPB;
        BNSN_P = BNSN_P / CNFLX;
        BNSN_C = BNSN_C / CNCCN;
        C1 = 2 * KUP * (KUP + 1);
        C2 = M - KUP - 1;
        %
        CADD = C1 / C2;
        AIC0(KUP) = CM * log(BNSNALL) + 2.0D0 * KUP + CADD; %! "c-AIC" by N. Sugiura (1978)
        AIC1(KUP) = CNAPB * log(BNSN_B) + CNFLX * log(BNSN_P) + 2.0D0 * KUP +CADD; %  ! "c-AIC" by N. Sugiura (1978)
        %
        % fprintf(% fid\d{3}, '%d %d %d %d %d %d %d\r\n', KUP, BNSNALL, BNSN_B, BNSN_P, ...
            % BNSN_C, AIC0(KUP), AIC1(KUP));
        % fprintf(% fid\d{3}, '%d %d\r\n', KUP, CNAPB * log(BNSN_B));
        % fprintf(% fid\d{3}, '%d %d\r\n', KUP, CNFLX * log(BNSN_P));
        % fprintf(% fid\d{3}, '%d %d\r\n', KUP, CNCCN * log(BNSN_C));
        % fprintf(% fid\d{3}, '%d %d\r\n', KUP, 2.0D0 * KUP + CADD);
        %
        if (AIC0(KUP) < AI0MIN)
            AI0MIN = AIC0(KUP);
            K0UP = KUP;
        end

        if (AIC1(KUP) < AI1MIN)
            AI1MIN = AIC1(KUP);
            K1UP = KUP;
        end

        %
    end

    % fprintf('%s %d\r\n', 'Min.AIC1 occurs at KUP=', K0UP);
    % fprintf('%s %d\r\n', 'Min.AIC1 occurs at KUP=', K1UP);
    % fprintf('%s\r\n', 'AICを過信せず、特異値の並びに');
    % fprintf('%s\r\n', '注目して打ち切り箇所を決定せよ');
    AIC0(1:N) = AIC0(1:N) - AI0MIN;
    AIC1(1:N) = AIC1(1:N) - AI1MIN;

    for KUP = 1:N
        % fprintf(% fid\d{3}, '%d %d\r\n', KUP, AIC0(KUP));
        % fprintf(% fid\d{3}, '%d %d\r\n', KUP, AIC1(KUP));
    end

    
end
