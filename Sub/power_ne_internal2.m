function [X_ne, F] = power_ne_internal2(J0_cell, Xi_mat, Phi_mat, V, vec_x0, ...
                                        vec_p, inode_mask)
% VERSION 2 of the power allocaiton algorithm for partially coordinated
% game. We also consider there are internal neighbor nodes. 
% 
% - inode_mask
% indicates whether the nodes is internal node. 
%
% - vec_p 
% contains the power allocation of the internal neighbor node.
% 
% [X_ne, F] = power_ne_internal2(J0_cell, Xi_mat, Phi_mat, V_vec, Xne0)
% 
%   where F is the objective function value.
%
% Find the power allocaiton via NE computed using the internal node
% coordination method.
%
% NOTE: User 1 is the internal node.

K = size(Phi_mat, 1);
Kle = K - sum(inode_mask) - 1;  % Number of leaf neighbor nodes

if nargin < 5
    % vec_x0 = ones(1, K) / K;
    vec_x0 = rand(1, K);
    vec_x0 = vec_x0 / sum(vec_x0);
end

% Construct the parameter $Delta$
Delta = zeros(K);
for k = 1:K
    for j = 1:K
        if j == k
            continue;
        end
        phi = Phi_mat(k, j);
        u_kj = [cos(phi), sin(phi)].';
        J = J0_cell{j};
        % invJ = inv(J);
        delta_kj = u_kj' * (J \ u_kj);
        Delta(k, j) = delta_kj;
    end
end

% The centralized minimization problem
my_fun = @(vec_x) func_f(vec_x, 1, J0_cell, Xi_mat, Phi_mat, Delta, V, vec_p, inode_mask);
my_fun_A = ones(1, K);
my_fun_b = 1;
my_fun_lb = zeros(K, 1);
my_opt = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
% my_opt0 = optimoptions(@fmincon,'Algorithm','interior-point');

vec_xk = fmincon(my_fun, vec_x0, my_fun_A, my_fun_b, [], [], my_fun_lb, [], [], my_opt);

if size(vec_xk, 1) > size(vec_xk, 2)
    vec_xk = vec_xk.';
end

% Construct the full power allocation
X_ne = zeros(K);
X_ne(1, :) = vec_xk;
for j = 2:Kle + 1
    % BEST RESPONSE T_jk of x_kj (for k === 1)
    x_kj = vec_xk(j);
    
    xi_kj = Xi_mat(1, j);
    
    phi = Phi_mat(j, 1);
    u_jk = [cos(phi) sin(phi)].';
    D_jk = u_jk * u_jk.';

    delta_jk = Delta(j, 1);
    
    J0_j = J0_cell{j};
    
    a1 = sqrt(trace(J0_j * J0_j - J0_j * D_jk * J0_j) * xi_kj) - det(J0_j) * sqrt(V);
    b1 = (det(J0_j) * delta_jk + trace((eye(2) - D_jk)*J0_j)) * xi_kj * sqrt(V);
    c1 = det(J0_j) * sqrt(V);

    x_jk = max(0, min(a1 * x_kj / (c1 + b1 * x_kj), 1));
    X_ne(j, 1) = x_jk;
    
end

F = zeros(1, K);
for k = 1:K
    F(k) = func_f(X_ne(k, :), k, J0_cell, Xi_mat, Phi_mat, Delta, V, vec_p, inode_mask);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = func_f(vec_xk, ...     % vector, power alloc of user k to others
                    k, ...          % The k-th agent
                    JA, ...         % 1 * Nms cell, J_e^A, EFIM of user k
                    Xi_ms2ms, ...   % Nms * Nms matrix, the channel
                    Phi_ms2ms, ...  % Nms * Nms matrix, the angle
                    Delta, ...      % Nms * Nms matrix, \delta_k,l
                    V, ...          % The trade-off parameter
                    vec_p, ...      % Power allocation of the internal neighbor nodes
                    inode_mask)     % Indicator of the internal neighbor nodes
%

% NOTE: The sum power and peak power constraints are assumed to be 1.
% (normalized)

Nms = size(Phi_ms2ms, 1);

J = JA{k};
for j = 1:Nms
    if j == k 
        continue;
    end
    
    phi = Phi_ms2ms(k, j);
    d_kj = [cos(phi) sin(phi)].';
    D_kj = d_kj * d_kj';
    delta_kj = Delta(k, j);
    xi_kj = Xi_ms2ms(k, j);
        
    if inode_mask(j) == 0
        % BEST RESPONSE T_jk of x_kj
        x_kj = vec_xk(j);

        u_jk = [cos(phi + pi) sin(phi + pi)].';
        D_jk = u_jk * u_jk.';

        delta_jk = Delta(j, k);

        J0_j = JA{j};

        a1 = sqrt(trace(J0_j * J0_j - J0_j * D_jk * J0_j) * xi_kj) - det(J0_j) * sqrt(V);
        b1 = (det(J0_j) * delta_jk + trace((eye(2) - D_jk)*J0_j)) * xi_kj * sqrt(V);
        c1 = det(J0_j) * sqrt(V);

        x_jk = max(0, min(a1 * x_kj / (c1 + b1 * x_kj), 1));
        % --
    else
        x_jk = vec_p(j);
    end
    
    % Two-way TOA scheme
    eff_Gamma_kj = vec_xk(j) * x_jk * xi_kj ...    % effecive SNR (gamma) for the bound
            / (1e-10 + vec_xk(j) + x_jk + vec_xk(j) * x_jk * xi_kj * delta_kj);
             % 1e-10 is to avoid the 0/0 situation
    
    J = J + eff_Gamma_kj * D_kj;
end

f = trace_inv(J) + V * sum(vec_xk);

end