function [Xne, Fspeb, Fcost] = power_game1(J0_cell, Xi_mat, Phi_mat, V_vec, Conn)
% PHASE TWO RESEARCH (LARGE NETWORK GAME)
% THIS IS MODIFIED FROM power_ne_np.m, UPDATED WITH Graph_connection (I.E.,
% PARTIALLY CONNECTED GRAPH)
% 
% [Xne, Fspeb, Fcost] = power_ne_np(J0_cell, Xi_mat, Phi_mat, V_vec, Xne0)
%
% OUTPUT PARAMETERS
%   Xne     Nms * Nms matrix
%   Fspeb   1 * Nms row vector
%   Fcost   1 * Nms row vector, the objective values (cost)


N_max = 8;
eps = 1e-9;

% UPDATED
Ntrack = 5;
X_array = zeros(Ntrack, N_max);
track_id = zeros(Ntrack, 1);
for i = 1:Ntrack
    track_id(i) = find(Conn(i, :) > 0, 1);
end
% --


Nms = size(Xi_mat, 1); % # of players (users)

% Preparation
Delta = zeros(Nms);
for k = 1:Nms
    for j = 1 : Nms
        if k == j
            continue;
        end
        phi = Phi_mat(k, j);
        Delta(k, j) = [cos(phi) sin(phi)] * inv(J0_cell{j}) * [cos(phi) sin(phi)].';
    end
end

Xne0 = ones(Nms) / (Nms - 1);    
Xne = Xne0;
Xne0 = Xne0 + 1;

my_opt = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');

Fspeb = zeros(1, Nms);
Fcost = zeros(1, Nms);

i = 0;
while norm(Xne0 - Xne, 'fro') > eps && i < N_max
    i = i + 1;
    
    Xne0 = Xne;
    for k = 1:Nms
        vec_x0 = Xne0(k, :).';
        vec_x1 = Xne0(:, k);
        
        my_fun = @(x) func_ispeb(x, vec_x1, k, J0_cell, Xi_mat, Phi_mat, ...
                                 Delta, Conn, V_vec(k));
        my_fun_A = ones(1, Nms);
        my_fun_b = 1;
        my_fun_lb = zeros(Nms, 1);
        [x_k, f_k] = fmincon(my_fun, vec_x0, my_fun_A, my_fun_b, [], [], my_fun_lb, [], [], my_opt);
        Xne(k, :) = x_k;
        Fcost(k) = f_k;
        Fspeb(k) = f_k - V_vec(k) * sum(x_k);
        
    end
    
    for j = 1:Ntrack
        X_array(j, i) = Xne(j, track_id(j));
    end
end

num_iter = i;
% if i < N_max
%     fprintf('NE found in %d steps.\n', i);
% else
%     fprintf('Not converged.\n');
% end

% % UPDATED
% figure(2),
% PowerBudget = 2; %<- outside parameter for consistency
% my_markers = {'d', 's', '^', 'o', '<', '>', 'v', '+', 'p', '.', '*'}.';
% curve_str = cell(Ntrack, 1);
% for i = 1:Ntrack
%     curve_str{i} = sprintf('Agent %d, x%d%d', i, i, track_id(i));
% end
% p_handle = plot(0:N_max, [ones(Ntrack, 1) / (Nms - 1) X_array] * PowerBudget, 'LineWidth', 2);
% set(p_handle, {'marker'}, my_markers(1:Ntrack), 'MarkerSize', 9);
% set(gca, 'FontSize', 14);
% legend(curve_str);
% ylabel('Power allocation');
% xlabel('Iteration');
% xlim([0, num_iter]);
% set(gca, 'XTick', 0:1:num_iter);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = func_ispeb(vec_x0, ...     % vector, power alloc of user k to others
                        vec_x1, ...     % vector, power alloc of others to user k
                        k, ...          % The k-th agent
                        JA, ...         % 1 * Nms cell, J_e^A, EFIM of user k
                        Xi_ms2ms, ...   % Nms * Nms matrix, the channel
                        Phi_ms2ms, ...  % Nms * Nms matrix, the angle
                        Delta, ...      % Nms * Nms matrix, \delta_k,l
                        Conn, ...       % [UPDATED] Graph connection
                        V)              % The trade-off parameter
%
Nms = size(Phi_ms2ms, 1);

J = JA{k};
for j = 1:Nms
    if j == k 
        continue;
    end
    if Conn(k, j) == 0  % UPDATED
        continue;
    end
    
    phi = Phi_ms2ms(k, j);
    d_kj = [cos(phi) sin(phi)].';
    D_kj = d_kj * d_kj';
    delta_kj = Delta(k, j);
    xi_kj = 4 * Xi_ms2ms(k, j);  % CAUTION: factor of 4!
    
%     % One-way TOA scheme
%     eff_gamma_kj = vec_x0(j) * xi_kj ...    % effecive SNR (gamma) for the bound
%             / (1 + vec_x0(j) * xi_kj * delta_kj);
        
    % Two-way TOA scheme
    eff_gamma_kj = vec_x0(j) * vec_x1(j) * xi_kj ...    % effecive SNR (gamma) for the bound
            / (1e-11 + vec_x0(j) + vec_x1(j) + vec_x0(j) * vec_x1(j) * xi_kj * delta_kj);
             % 1e-11 is to avoid the 0/0 situation
    
    J = J + eff_gamma_kj * D_kj;
end

f = trace_inv(J) + V * sum(vec_x0);

end