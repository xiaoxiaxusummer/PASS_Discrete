function F = optimize_F_HybridMIMO(K, N, M, W, F, D, lambda_W, rho, subconnected)

C = zeros(N,N)*1j; 
B = zeros(M,N)*1j;

for k = 1:K
    C = C + D(:,k) * (D(:,k)');
    B = B + (W(:,k) + rho * lambda_W(:,k)) * (D(:,k)');
end

A = eye(M);
F = opt_BCD(F, A, C, B, rho, W, D, lambda_W, subconnected);

end


function X = opt_BCD(X_init, A, C, B, rho, W, D, lambda_W, subconnected)
    tol = 1e-6;
    max_iter = 100;
    [M, N] = size(X_init);
    L = round(M/N);
    X = X_init;
    Z = A * X * C;
    for iter = 1:max_iter
        X_prev = X;
        for i = 1:M
            for j = 1:N
                if subconnected && ((i<(j-1)*L+1) || (i>j*L))
                    X(i,j) = 0;
                else
                    b = A(i,i) *X(i,j) *C(j,j) - Z(i,j) + B(i,j);
                    v = b/abs(b);
                    Z = Z + (v - X(i,j)) * (A(:,i) * (C(j,:)));
                    X(i,j) = v;
                end
            end
        end
        residual_W = W - X*D;
        % disp(1/2/rho*sum(sum_square_abs(residual_W)));
        L_W = 1/2/rho * sum(sum(abs(residual_W + lambda_W * rho).^2));
        if sum(sum_square_abs(X-X_prev)) < tol
            break;
        end
    end
end