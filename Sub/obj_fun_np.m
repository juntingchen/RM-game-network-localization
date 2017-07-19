function [Fspeb, Fcost] = obj_fun_np(X, J0_cell, Xi_mat, Phi_mat, V)
% VERSION 2: fixed the channel gain issue by factor 4

Xi_mat = Xi_mat * 4;

K = size(Xi_mat, 1);

Fspeb = zeros(1, K);
Fcost = zeros(1, K);
for k = 1:K
    J = J0_cell{k};
    
    for j = 1:K
        if j == k 
            continue;
        end

        phi = Phi_mat(k, j);
        u_kj = [cos(phi) sin(phi)].';
        D_kj = u_kj * u_kj';
        delta_kj = u_kj' * inv(J0_cell{j}) * u_kj;
        xi_kj = Xi_mat(k, j);

        eff_gamma_kj = X(k,j) * X(j,k) * xi_kj ...    % effecive SNR (gamma) for the bound
                / (1e-20 + X(k, j) + X(j, k) + X(k, j) * X(j, k) * xi_kj * delta_kj);
                 % 1e-11 is to avoid the 0/0 situation

        J = J + eff_gamma_kj * D_kj;
    end

    Fspeb(k) = trace(inv(J));
    Fcost(k) = Fspeb(k) + V * sum(X(k, :));
end