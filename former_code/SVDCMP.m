% function [W,V,A] = SVDCMP(A,M,N,MP,NP)
function [W, V, A] = SVDCMP(A)
    %%%NMAX = 100;
    NMAX = 4900; %    !12月5日 NNEQと連動
    [M, N] = size(A);
    %A = zeros(MP,NP);
    % W = zeros(1,NP);
    % V = zeros(NP,NP);
    W = zeros(1, N);
    V = zeros(N, N);
    RV1 = zeros(1, NMAX);
    %
    if (NMAX >= N)

        if (M < N)
            fprintf('%s\r\n', 'You must augment A with extra zero rows.');
            pause
        end

        %    ! Householder reduction to bidiagonal form.
        G = 0.0;
        SCALE = 0.0;
        ANORM = 0.0;

        for I = 1:N
            L = I + 1;
            RV1(I) = SCALE * G;
            G = 0.0;
            S = 0.0;
            SCALE = 0.0;

            if (I <= M)

                for K = 1:M
                    SCALE = SCALE + abs(A(K, I));
                end

                if (SCALE ~= 0.0)

                    for K = I:M
                        A(K, I) = A(K, I) / SCALE;
                        S = S + A(K, I) * A(K, I);
                    end

                    F = A(I, I);
                    G = -abs(sqrt(S)) * sign(F); %-DSIGN(sqrt(S),F)
                    H = F * G - S;
                    A(I, I) = F - G;

                    if (I ~= N)

                        for J = L:N
                            S = 0.0;

                            for K = I:M
                                S = S + A(K, I) * A(K, J);
                            end

                            F = S / H;

                            for K = I:M
                                A(K, J) = A(K, J) + F * A(K, I);
                            end

                        end

                    end

                    for K = I:M
                        A(K, I) = SCALE * A(K, I);
                    end

                end

            end

            W(I) = SCALE * G;
            G = 0.0;
            S = 0.0;
            SCALE = 0.0;

            if and(I <= M, I ~= N)

                for K = L:N
                    SCALE = SCALE + abs(A(I, K));
                end

                if (SCALE ~= 0.0)

                    for K = L:N
                        A(I, K) = A(I, K) / SCALE;
                        S = S + A(I, K) * A(I, K);
                    end

                    F = A(I, L);
                    G = -abs(sqrt(S)) * sign(F);
                    H = F * G - S;
                    A(I, L) = F - G;

                    for K = L:N
                        RV1(K) = A(I, K) / H;
                    end

                    if (I ~= M)

                        for J = L:M
                            S = 0.0;

                            for K = L:N
                                S = S + A(J, K) * A(I, K);
                            end

                            for K = L:N
                                A(J, K) = A(J, K) + S * RV1(K);
                            end

                        end

                    end

                    for K = L:N
                        A(I, K) = SCALE * A(I, K);
                    end

                end

            end

            ANORM = max(ANORM, (abs(W(I)) + abs(RV1(I))));
        end

        %! Accumulation of right-hand transformations.
        for I = N:-1:1

            if (I < N)

                if (G ~= 0.0)

                    for J = L:N
                        V(J, I) = (A(I, J) / A(I, L)) / G;
                    end

                    for J = L:N
                        S = 0.0;

                        for K = L:N
                            S = S + A(I, K) * V(K, J);
                        end

                        for K = L:N
                            V(K, J) = V(K, J) + S * V(K, I);
                        end

                    end

                end

                for J = L:N
                    V(I, J) = 0.0;
                    V(J, I) = 0.0;
                end

            end

            V(I, I) = 1.0;
            G = RV1(I);
            L = I;
        end

        %    ! Accumulation of left-hand transformations.
        for I = N:-1:1
            L = I + 1;
            G = W(I);

            if (I < N)

                for J = L:N
                    A(I, J) = 0.0;
                end

            end

            if (G ~= 0.0)
                G = 1.0 / G;

                if (I ~= N)

                    for J = L:N
                        S = 0.0;

                        for K = L:M
                            S = S + A(K, I) * A(K, J);
                        end

                        F = (S / A(I, I)) * G;

                        for K = I:M
                            A(K, J) = A(K, J) + F * A(K, I);
                        end

                    end

                end

                for J = I:M
                    A(J, I) = A(J, I) * G;
                end

            else

                for J = I:M
                    A(J, I) = 0.0;
                end

            end

            A(I, I) = A(I, I) + 1.0;
        end

        %   ! Diagonalization of the bidiagonal form.
        for K = N:-1:1

            for ITS = 1:300

                for L = K:-1:1
                    NM = L - 1;

                    if ((abs(RV1(L)) + ANORM) == ANORM)
                        break
                    elseif ((abs(W(NM)) + ANORM) == ANORM)
                        C = 0.0;
                        S = 1.0;

                        for I = L:K
                            F = S * RV1(I);

                            if ((abs(F) + ANORM) ~= ANORM)
                                G = W(I);
                                H = sqrt(F * F + G * G);
                                W(I) = H;
                                H = 1.0 / H;
                                C = (G * H);
                                S = -(F * H);

                                for J = 1:M
                                    Y = A(J, NM);
                                    Z = A(J, I);
                                    A(J, NM) = (Y * C) + (Z * S);
                                    A(J, I) = -(Y * S) + (Z * C);
                                end

                            end

                        end

                        break
                    end

                end

                Z = W(K);

                if (L == K)

                    if (Z < 0.0)
                        W(K) = -Z;

                        for J = 1:N
                            V(J, K) = -V(J, K);
                        end

                    end

                    break
                end

                if (ITS == 300)
                    fprintf('%s\r\n', 'No convergence in 300 iterations');
                    pause
                end

                X = W(L);
                NM = K - 1;
                Y = W(NM);
                G = RV1(NM);
                H = RV1(K);
                F = ((Y - Z) * (Y + Z) + (G - H) * (G + H)) / (2.0 * H * Y);
                G = sqrt(F * F + 1.0);
                F = ((X - Z) * (X + Z) + H * ((Y / (F + abs(G) * sign(F))) - H)) / X;
                %          ! Next QR transformation:
                C = 1.0;
                S = 1.0;

                for J = L:NM
                    I = J + 1;
                    G = RV1(I);
                    Y = W(I);
                    H = S * G;
                    G = C * G;
                    Z = sqrt(F * F + H * H);
                    RV1(J) = Z;
                    C = F / Z;
                    S = H / Z;
                    F = (X * C) + (G * S);
                    G = -(X * S) + (G * C);
                    H = Y * S;
                    Y = Y * C;

                    for JJ = 1:N
                        X = V(JJ, J);
                        Z = V(JJ, I);
                        V(JJ, J) = (X * C) + (Z * S);
                        V(JJ, I) = -(X * S) + (Z * C);
                    end

                    Z = sqrt(F * F + H * H);
                    W(J) = Z;
                    %                if (Z ~= 0.0)
                    Z = (1.0 / Z) .* (Z ~= 0.0) + Z .* (Z == 0.0);
                    C = F * Z .* (Z ~= 0.0) + C .* (Z == 0.0);
                    S = H * Z .* (Z ~= 0.0) + S .* (Z == 0.0);
                    %                end
                    F = (C * G) + (S * Y);
                    X = -(S * G) + (C * Y);

                    for JJ = 1:M
                        Y = A(JJ, J);
                        Z = A(JJ, I);
                        A(JJ, J) = (Y * C) + (Z * S);
                        A(JJ, I) = -(Y * S) + (Z * C);
                    end

                end

                RV1(L) = 0.0;
                RV1(K) = F;
                W(K) = X;
            end

            %    3 CONTINUE
        end

    else
        fprintf('%s\r\n', 'Lack of Dimension (NMAX<N) in Subr. SVBKSB');
    end
end
