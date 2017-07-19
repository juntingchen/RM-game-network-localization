function [Xne, Fspeb, Fcost] = power_game4(J0_cell, Xi_ms2ms, Phi_ms2ms, ...
                                            V0_vec, Pt_vec, Conn)
% [Xne1, Fcost] = power_game3(J0_cell, Xi_ms2ms, Phi_ms2ms, V_vec, Xne0)
%
% THIS IS BUILT FROM POWER_GAME3.m 
%
% OUTPUT PARAMETERS
%   Xne1    Nms * Nms matrix
%   Fspeb   1 * Nms row vector
%   Fcost   1 * Nms row vector, the objective values (cost)


% Ntrack = 5; % <- for debug purpose

Veps = mean(V0_vec(:)) * 1e-2;

Nms = size(Xi_ms2ms, 1); % # of players (users)
% V = V0_vec(1);

% Preparation
Delta = zeros(Nms);
Cmat = ones(Nms);
for k = 1:Nms
    for j = 1 : Nms
        if k == j
            continue;
        end
        phi = Phi_ms2ms(k, j);
        Delta(k, j) = [cos(phi) sin(phi)] * inv(J0_cell{j}) * [cos(phi) sin(phi)].';
        Cmat(k, j) =[cos(phi) sin(phi)] * (inv(J0_cell{k})^2) * [cos(phi) sin(phi)].';
    end
end
Delta2 = Delta + Delta.';

% --------
MAXLOOP = 20;

% Vs_array = zeros(Ntrack, MAXLOOP);  % <- for debug purpose

Xne = zeros(Nms);
vt_vec = V0_vec;
i_out = 0;
vt_gap = 1;
while i_out < MAXLOOP && vt_gap > Veps
    i_out = i_out + 1;
   
    Gamma = repmat(vt_vec(:), 1, Nms) ./ Cmat;
    
    vt1_vec = zeros(1, Nms);
    for k = 1:Nms
        Peps = Pt_vec(k) * 3e-2;
        Gamma_k = Gamma;
        
        Vk = vt_vec(k);
        Vk_min = V0_vec(k);
        Vk_max = Vk;
        
        i_in = 0;
        Pksum = 0;
        Vk0 = - 1;
        while i_in < MAXLOOP && abs(Pksum - Pt_vec(k)) > Peps ...
                && ~(Pksum < Pt_vec(k) && abs(Vk0 - V0_vec(k)) < Veps / 2)
            i_in = i_in + 1;
            Xk = zeros(1, Nms);
            
            Vk0 = Vk;
            Gamma_k(k, :) = Vk * ones(1, Nms) ./ Cmat(k, :);
            for j = 1:Nms
                if Conn(k, j) == 0
                    continue;
                end
                akj = 2 * sqrt(Xi_ms2ms(k, j) / Gamma_k(k, j)) - 1;
                ajk = 2 * sqrt(Xi_ms2ms(k, j) / Gamma_k(j, k)) - 1;
                bkj = 4 * Xi_ms2ms(k, j) * Delta2(k, j);
                Xk(j) = max(0, (akj * ajk - 1) / (bkj * (1 + ajk)));
            end
            Pksum = sum(Xk);
            
            if Pksum < Pt_vec(k) - Peps
                Vk_max = Vk;
                Vk1 = (Vk + Vk_min) / 2;
                Vk = Vk1;
                % if abs(Vk1 - V0_vec(k)) < Veps
            elseif Pksum > Pt_vec(k) + Peps
                if abs(Vk - Vk_max) < Veps / 2
                    Vk_max = Vk_max * 2;
                    Vk_min = Vk;
                    Vk1 = Vk * 2;
                else
                    Vk_min = Vk;
                    Vk1 = (Vk + Vk_max) / 2;
                end
                Vk = Vk1;
            end
            
                    
        end
        Xne(k, :) = Xk;
        vt1_vec(k) = Vk;
        Gamma(k, :) = Gamma_k(k, :);
        
    end
    vt_gap = max(abs(vt1_vec - vt_vec));
    vt_vec = vt1_vec;
    
    % Vs_array(:, i_out) = vt_vec(1 : Ntrack);

end

% % For convergence test
% figure(1),
% my_markers = {'d', 's', '^', 'o', '<', '>', 'v', '+', 'p', '.', '*'}.';
% curve_str = cell(Ntrack, 1);
% for i = 1:Ntrack
%     curve_str{i} = sprintf('Agent %d', i);
% end
% p_handle = plot(0:MAXLOOP, [V0_vec(1:Ntrack).' Vs_array], 'LineWidth', 2);
% set(p_handle, {'marker'}, my_markers(1:Ntrack), 'MarkerSize', 9);
% set(gca, 'FontSize', 14);
% legend(curve_str);
% ylabel('V');
% xlabel('Iteration');
% xlim([0, i_out]);
% set(gca, 'XTick', 0:1:i_out);
% return;


Fspeb = zeros(1, Nms);
Fcost = zeros(1, Nms);

for k = 1:Nms
    vec_x0 = Xne(k, :).';
    vec_x1 = Xne(:, k);
    Fcost(k) = func_f(vec_x0, vec_x1, k, J0_cell, Xi_ms2ms, Phi_ms2ms, Delta, V0_vec(k));
    Fspeb(k) = Fcost(k) - V0_vec(k) * sum(Xne(k, :));
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = func_f(vec_x0, ...     % vector, power alloc of user k to others
                    vec_x1, ...     % vector, power alloc of others to user k
                    k, ...          % The k-th agent
                    JA, ...         % 1 * Nms cell, J_e^A, EFIM of user k
                    Xi_ms2ms, ...   % Nms * Nms matrix, the channel
                    Phi_ms2ms, ...  % Nms * Nms matrix, the angle
                    Delta, ...      % Nms * Nms matrix, \delta_k,l
                    V)              % The trade-off parameter
%
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

    %     % One-way TOA scheme
    %     eff_gamma_kj = vec_x0(j) * xi_kj ...    % effecive SNR (gamma) for the bound
    %             / (1 + vec_x0(j) * xi_kj * delta_kj);

        % Two-way TOA scheme
        eff_gamma_kj = vec_x0(j) * vec_x1(j) * xi_kj ...    % effecive SNR (gamma) for the bound
                / (1e-11 + vec_x0(j) + vec_x1(j) + vec_x0(j) * vec_x1(j) * xi_kj * delta_kj);
                 % 1e-11 is to avoid the 0/0 situation

        J = J + eff_gamma_kj * D_kj;
    end

    f = trace(inv(J)) + V * sum(vec_x0);

end

% ----
function F = func_nbs(Fd,...          % 
                      vec_X, ...      % 
                      J0_cell, ...    % 1 * Nms cell, J_e^A
                      Xi_ms2ms, ...   % Nms * Nms matrix, the channel
                      Phi_ms2ms, ...  % Nms * Nms matrix, the angle
                      Delta, ...      % Nms * Nms matrix, \delta_k,l
                      V,...           % The trade-off parameter
                      f_mask)         % Feasibility mask
%
    Nms = size(Phi_ms2ms, 1);
    X = reshape(vec_X, Nms, Nms);
    
    F = 0;
    for k = 1:Nms
        if f_mask(k) > 0
            vec_x0 = X(k, :).';
            vec_x1 = X(:, k);
            fk = func_f(vec_x0, vec_x1, k, J0_cell, Xi_ms2ms, Phi_ms2ms, Delta, V);
            F = F - log(Fd(k) - fk);
        end
    end

end
% ------
% Constraint function for the Feasibility problem
function [C, ceq] = myconfun_t(F0, vec_X, J0_cell, Xi_ms2ms, Phi_ms2ms, Delta, V)
    Nms = size(Xi_ms2ms, 1);
    X = reshape(vec_X(1:end-1), Nms, Nms);
    t = vec_X(end);
    
    C = zeros(Nms, 1);
    for k = 1:Nms
        vec_x0 = X(k, :).';
        vec_x1 = X(:, k);
        C(k) = func_f(vec_x0, vec_x1, k, J0_cell, Xi_ms2ms, Phi_ms2ms, Delta, V)...
               - F0(k) - t;
    end
    ceq = 0;
end

% ------
% Constraint function for the NBS problem
function [C, ceq] = myconfun(F0, vec_X, J0_cell, Xi_ms2ms, Phi_ms2ms, Delta, V)
    Nms = size(Xi_ms2ms, 1);
    X = reshape(vec_X, Nms, Nms);
    
    C = zeros(Nms, 1);
    for k = 1:Nms
        vec_x0 = X(k, :).';
        vec_x1 = X(:, k);
        C(k) = func_f(vec_x0, vec_x1, k, J0_cell, Xi_ms2ms, Phi_ms2ms, Delta, V)...
               - F0(k);
    end
    ceq = 0;
end