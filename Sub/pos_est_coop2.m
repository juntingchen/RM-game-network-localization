function p_hat_g1 = pos_est_coop2(Xne1, p_hat, J0_cell, Dist_ms2ms, Xi_ms2ms, Phi_ms2ms)
%
% Version 2: using prior informaiton, i.e., ls_coop -> ls_coop_map
%
% Note 1: this version has fixed the chanel gain issue by a factor of 4
% for asynchronous networks
% 
% Note 2: Pt_bs * Xi_ms2bs --> Xi_ms2bs

OPTIONS = optimset('Algorithm','levenberg-marquardt', ...
                           'Display', 'off');
                       
Nms = size(Xi_ms2ms, 1);

p_hat_g1 = zeros(Nms, 2);   % Position estimation
for i_ms = 1:Nms 
    cnt = 0;   
    Z_c = zeros(1, Nms - 1);
    Lambda_c = zeros(1, Nms - 1);
    Pos_hat = zeros(Nms - 1, 2);
    for j = 1:Nms
        if j == i_ms
            continue;
        end
        cnt = cnt + 1;
        gamma_c = 4 * Xne1(i_ms, j) * Xne1(j, i_ms) * Xi_ms2ms(i_ms, j) ...
                       / (1e-10 + Xne1(i_ms, j) + Xne1(j, i_ms));
                   % Note that the factor 4 is to fix the channel gain
                   % issue for asynchronous networks
        Z_c(cnt) = Dist_ms2ms(i_ms, j) + randn * sqrt(1 / (1e-8 + gamma_c));
                                        % Inter-agent ranging measurements
        delta = [cos(Phi_ms2ms(i_ms, j)) sin(Phi_ms2ms(i_ms, j))] ...
                * inv(J0_cell{j}) * [cos(Phi_ms2ms(i_ms, j)) sin(Phi_ms2ms(i_ms, j))]';
        Lambda_c(cnt) = gamma_c / (1 + gamma_c * delta);
        Pos_hat(cnt, :) = p_hat(j, :);
    end

    if norm(Lambda_c) > 0
        ls_func = @(x) ls_coop_map(x, p_hat(i_ms, :), J0_cell{i_ms},  ...
                               Lambda_c, Z_c, Pos_hat);
        p_hat_g1(i_ms, :) = lsqnonlin(ls_func, p_hat(i_ms, :), [], [], OPTIONS);
    else
        p_hat_g1(i_ms, :) = p_hat(i_ms, :);
    end
end