% This script generates simulation results depicted in Fig. 3 of the paper
%
% J. Chen, W. Dai, Y. Shen, and M. Z. Win, "Resource Management Games for 
% Distributed Network Localization," IEEE J. Sel. Areas Commun., vol. 35, 
% no. 2, pp. 317--329, Feb. 2017.
%
% Author: Junting Chen (juntingc@usc.edu)
% Date: Sept 19, 2017.
%
% The goal of this simulation is to evaluate the SPEB performance and the
% power consumption over network cooperative localization with round-trip
% TOA ranging under various game approaches.
%
% ABBREVIATIONS
%   ERC     Equivalent ranging coefficient
%   FIM     Fisher information matrix
%
% ANNOTATION
%   Game 4 stands for the per-link negotiation game in the paper.
%
% AUGMENTED Simulation flow:
%   Generate topology (anchor positions)
%   (To repeat)
%   Generate channel quality (anchors, agents)
%   Compute ERC and the initial FIM
%   Compute power allocation for cooperation
%   Estimation 1: Based on anchors only
%   Estimation 2: Cooperation
%   Performance evaluation
%   (Repeat end)
clear 
addpath Sub

% Case Configuration
Vs = 0.02;         % Changed from an array to a scaler
Sigma_dB = 3:3:30; % dBm 
N_topo = 2000;     % # of realizations (Monte-Carlo simulations)
L0 = 2;

% Power and Channel SETTING
Pt = 1; 
Pt_bs_dBm = 33;     % 2W
Pt_bs = 10 ^ (Pt_bs_dBm / 10);

Kfactor_dB = 4.7;
Kfactor = 10^(Kfactor_dB / 10);
fc = 5.25e9;		% Carrier frequency (pathloss model)

los_ratio = 0.05;
rice_nu = sqrt(Kfactor / (1 + Kfactor));
rice_s = sqrt(1 / 2 / (1 + Kfactor));

% Topology Setting
Random_Agent_Pos = 0;
NumberOfAgents = 36;     % This parameter will be suppressed when Random_Agent_Pos = 0
AreaLength = 200;

Scheme_str = {'No cooperation'
              'SE, L0 = 1'
              sprintf('LBE, L0 = %d', L0)
              'Centralized NE'
              'Naive scheme'
              };
          
% ----------------------------------------------------------------------- %
% ------------------------ Anchor Topology ------------------------------ %
% Placement of the anchors
Pos_bs = [1  -1
         -1   1
         -1  -1
          1   1
          ] * AreaLength / 6;
Nbs = size(Pos_bs, 1);
% ----

Sigma = 10.^(Sigma_dB / 10);
N_sigma = length(Sigma);
N_Vs = length(Vs);              % = 1 for this script

Speb0 = zeros(N_topo, N_sigma);      % SPEB
Speb1 = zeros(N_topo, N_sigma);
Speb_g1 = zeros(N_topo, N_sigma);
Speb_g2 = zeros(N_topo, N_sigma);
Speb_g3 = zeros(N_topo, N_sigma);
Speb_ms1_g1 = zeros(N_topo, N_sigma);
Speb_ms1_g2 = zeros(N_topo, N_sigma);
Speb_ms1_g3 = zeros(N_topo, N_sigma);
Speb_ms1_0 = zeros(N_topo, N_sigma);
Speb_ms1_1 = zeros(N_topo, N_sigma);

Pav_g1 = zeros(N_topo, N_sigma);	% Power consumption
Pav_g2 = zeros(N_topo, N_sigma);
Pav_g3 = zeros(N_topo, N_sigma);
Pav_ms1_g1 = zeros(N_topo, N_sigma);
Pav_ms1_g2 = zeros(N_topo, N_sigma);
Pav_ms1_g3 = zeros(N_topo, N_sigma);

Mse_nonc = zeros(N_topo, N_sigma);   % MSE, no cooperation
Mse_full = zeros(N_topo, N_sigma);	% MSE, cooperation, full power allocation
Mse_g1 = zeros(N_topo, N_sigma);   % MSE, cooperation
Mse_g2 = zeros(N_topo, N_sigma);   % MSE, cooperation
Mse_g3 = zeros(N_topo, N_sigma);   % MSE, cooperation
Mse_ms1_g1 = zeros(N_topo, N_sigma);   % MSE, cooperation
Mse_ms1_g2 = zeros(N_topo, N_sigma);   % MSE, cooperation
Mse_ms1_g3 = zeros(N_topo, N_sigma);   % MSE, cooperation
Mse_ms1_nonc = zeros(N_topo, N_sigma);
Mse_ms1_full = zeros(N_topo, N_sigma);

LinkNum = zeros(N_topo, 4);
parfor i_topo = 1:N_topo     % <- parfor HERE
    fprintf('%d\n', i_topo);
    % Eeah realization has an independnt channel and/or agent deployment
    
    % ---------------- Agent topology, channel, EFIM, etc. ------------------ %
    % Placement of the mobiles
    if ~Random_Agent_Pos
        Nms = NumberOfAgents;
        SquareRootNumberOfAgents = round(sqrt(NumberOfAgents));
        gridwidth = AreaLength / (SquareRootNumberOfAgents + 1);
        Pos_ms = zeros(Nms, 2);
        for i = 1:Nms
            row_id = floor((i - 1) / SquareRootNumberOfAgents) + 1;
            col_id = i - (row_id - 1) * SquareRootNumberOfAgents;
            
            posx = gridwidth * col_id - AreaLength / 2;
            posy = gridwidth * row_id - AreaLength / 2;
            Pos_ms(i, :) = [posx, posy];
        end
        % Pos_ms(randperm(Nms), :) = Pos_ms;
    else
        Nms = NumberOfAgents;
        Pos_ms = zeros(Nms, 2);
        Pos_ms(1, :) = [0 0];   % The 1st agent is set at origin
        for i = 2:Nms
            posx = rand * AreaLength - AreaLength / 2;
            posy = rand * AreaLength - AreaLength / 2;
            Pos_ms(i, :) = [posx, posy];
        end
    end
    if i_topo == 1
        hf1 = plot_topology(Pos_bs, Pos_ms, 1, AreaLength);
        axis square
    end
    
    V_vec = Vs * ones(1, Nms);

    % Parameters for EFIM
    Phi_ms2bs = zeros(Nms, Nbs);    
    Dist_ms2bs = zeros(Nms, Nbs);
    for i_ms = 1:Nms
        for i_bs = 1:Nbs
            vec = Pos_bs(i_bs, :) - Pos_ms(i_ms, :);
            Phi_ms2bs(i_ms, i_bs) = atan2(vec(2), vec(1));
            Dist_ms2bs(i_ms, i_bs) = norm(vec, 2);
        end
    end

    Phi_ms2ms = zeros(Nms);
    Dist_ms2ms = zeros(Nms);
    for k = 1:Nms
        for j = 1:Nms
            if k == j
                continue;
            end
            Dist_ms2ms(k, j) = norm(Pos_ms(k, :) - Pos_ms(j, :), 2);
            Dist_ms2ms(j, k) = Dist_ms2ms(k, j);
            % Phi_ms2ms(k, j) = atan((Pos_ms(j, 2) - Pos_ms(k, 2)) / (Pos_ms(j, 1) - Pos_ms(k, 1)));
            Phi_ms2ms(k, j) = atan2((Pos_ms(j, 2) - Pos_ms(k, 2)), (Pos_ms(j, 1) - Pos_ms(k, 1))); % <- suggested use
        end
    end
    
    [Xi0_ms2bs, Xi0_ms2ms, Snr0_ms2ms] = gen_chann(Dist_ms2bs, Dist_ms2ms);
    
%     mean(Xi0_ms2bs(:) * 1000),
%     mean(Xi0_ms2ms(1, :) * 100),
%     sort(Snr0_ms2ms(1, :)),
    
    Xi_ms2bs = Xi0_ms2bs;
    
    % Compute the initial FIM  
    J0_cell = cell(Nms, 1);
    for i_ms = 1:Nms
        J0_cell{i_ms} = zeros(2);
        for i_bs = 1:Nbs
            J0_cell{i_ms} = J0_cell{i_ms} + Xi_ms2bs(i_ms, i_bs) * Pt_bs ...
                * [cos(Phi_ms2bs(i_ms, i_bs)) sin(Phi_ms2bs(i_ms, i_bs))]' ...
                * [cos(Phi_ms2bs(i_ms, i_bs)) sin(Phi_ms2bs(i_ms, i_bs))];
        end
    end

    % Anchors' measurements
    OPTIONS = optimset('Algorithm','levenberg-marquardt', ...
                       'Display', 'off');
    p_hat = zeros(Nms, 2);
    Z_bs = Dist_ms2bs ...
            + randn(Nms, Nbs) .* sqrt(1 ./ (Xi_ms2bs * Pt_bs)); % measurement noise
    serr_non_temp = zeros(Nms, 1);
    for i_ms = 1:Nms
        ls_func = @(x) ls_noncoop(x, Pt_bs * Xi_ms2bs(i_ms, :), Z_bs(i_ms, :), Pos_bs);
        p_hat(i_ms, :) = lsqnonlin(ls_func, Pos_ms(i_ms,:), [], [], OPTIONS);
        serr_non_temp(i_ms) = norm(p_hat(i_ms, :) - Pos_ms(i_ms, :))^2;
    end
    Mse_nonc(i_topo, :) = mean(serr_non_temp) * ones(1, N_sigma);
    Mse_ms1_nonc(i_topo, :) = serr_non_temp(1) * ones(1, N_sigma);
    % -- end anchor's measurement --
    % ----------------------------------------------------------------------- %       
    
    numLinks = zeros(N_sigma, 4);
    for i_sigma = 1:N_sigma
        sigma = Sigma(i_sigma);
        Xi_ms2ms = Xi0_ms2ms * sigma;

        % GAME 1: Tree topology
        Xne1 = power_game2(J0_cell, Xi_ms2ms, Phi_ms2ms, V_vec);
        [Fspeb1, ~] = obj_fun_np(Xne1, J0_cell, Xi_ms2ms, Phi_ms2ms, Vs(1));
        Pav_g1(i_topo, i_sigma) = mean(sum(Xne1), 2);   % Average power consupmtion over all users
        Pav_ms1_g1(i_topo, i_sigma) = sum(Xne1(1, :));
        Speb_g1(i_topo, i_sigma) = mean(Fspeb1); 
        Speb_ms1_g1(i_topo, i_sigma)= Fspeb1(1);

        % cooperative position estimation
%         p_hat_g1 = pos_est_coop(Xne1, p_hat, J0_cell, Dist_ms2ms, Xi_ms2ms, Phi_ms2ms, ...
%                   Pt_bs * Xi_ms2bs, Z_bs, Pos_bs);
        p_hat_g1 = pos_est_coop2(Xne1, p_hat, J0_cell, Dist_ms2ms, Xi_ms2ms, Phi_ms2ms);
        serr = zeros(1, Nms);
        for i_ms = 1:Nms
            serr(i_ms) = norm(p_hat_g1(i_ms, :) - Pos_ms(i_ms, :))^2;
        end
        Mse_g1(i_topo, i_sigma) = mean(serr);
        Mse_ms1_g1(i_topo, i_sigma) = serr(1);
        
        numLinks(i_sigma, 1) = length(find(Xne1 > 0)) / 2;
        % -------- end game 1 --------

        % GAME 2: Per-link Negociation game (L = L0)
        Graph = pg_link_sel(J0_cell, Xi_ms2ms, Phi_ms2ms, V_vec, L0);
        Xne2 = power_game4(J0_cell, Xi_ms2ms, Phi_ms2ms, V_vec, Pt * ones(1, Nms), Graph);
        [Fspeb2, ~] = obj_fun_np(Xne2, J0_cell, Xi_ms2ms, Phi_ms2ms, Vs(1));
        Pav_g2(i_topo, i_sigma) = mean(sum(Xne2), 2);   % Average power consupmtion over all users
        Pav_ms1_g2(i_topo, i_sigma) = sum(Xne2(1, :));
        Speb_g2(i_topo, i_sigma) = mean(Fspeb2); 
        Speb_ms1_g2(i_topo, i_sigma) = Fspeb2(1);

%         p_hat_g2 = pos_est_coop(Xne2, p_hat, J0_cell, Dist_ms2ms, Xi_ms2ms, Phi_ms2ms, ...
%                   Pt_bs * Xi_ms2bs, Z_bs, Pos_bs);
        p_hat_g2 = pos_est_coop2(Xne2, p_hat, J0_cell, Dist_ms2ms, Xi_ms2ms, Phi_ms2ms);
        serr = zeros(1, Nms);
        for i_ms = 1:Nms
            serr(i_ms) = norm(p_hat_g2(i_ms, :) - Pos_ms(i_ms, :))^2;
        end
        Mse_g2(i_topo, i_sigma) = mean(serr);
        Mse_ms1_g2(i_topo, i_sigma) = serr(1);
        
        numLinks(i_sigma, 2) = length(find(Xne2 > 0)) / 2;
        % -- end game 2 --

        % GAME 3: Centralized NE
        [Graph1, NodeDegs, GpDegs] = pg_link_sel(J0_cell, Xi_ms2ms, Phi_ms2ms, V_vec, 8);
        Xne3 = power_game1(J0_cell, Xi_ms2ms, Phi_ms2ms, V_vec, Graph1);
            % This is for intilizing a good initial point for better
            % convergence (demonstration purpose)
        
        [Fspeb3, ~] = obj_fun_np(Xne3, J0_cell, Xi_ms2ms, Phi_ms2ms, Vs(1));
        Pav_g3(i_topo, i_sigma) = mean(sum(Xne3), 2);   % Average power consupmtion over all users
        Pav_ms1_g3(i_topo, i_sigma) = sum(Xne3(1, :));
        Speb_g3(i_topo, i_sigma) = mean(Fspeb3); 
        Speb_ms1_g3(i_topo, i_sigma) = Fspeb3(1);

        p_hat_g3 = pos_est_coop(Xne3, p_hat, J0_cell, Dist_ms2ms, Xi_ms2ms, Phi_ms2ms, ...
                  Pt_bs * Xi_ms2bs, Z_bs, Pos_bs);
        serr = zeros(1, Nms);
        for i_ms = 1:Nms
            serr(i_ms) = norm(p_hat_g3(i_ms, :) - Pos_ms(i_ms, :))^2;
        end
        Mse_g3(i_topo, i_sigma) = mean(serr);
        Mse_ms1_g3(i_topo, i_sigma) = serr(1);

          numLinks(i_sigma, 3) = length(find(Xne3 > 0)) / 2;
        % -- end game 3 --
        
		% Cooperative localization (NAIVE FULL power allocation)
        % Xfull = Pt / (Nms - 1) * (ones(Nms) - eye(Nms));
        Snr_ms2ms = sigma * Snr0_ms2ms;
        LinkSelectionThreshold_dB = 10;
        Sel = zeros(size(Snr_ms2ms));
        Sel(Snr_ms2ms > 10^(LinkSelectionThreshold_dB / 10)) = 1;
        Xfull = zeros(Nms);
        for i = 1:Nms
            
            cntLinks = length(find(Sel(i, :) > 0));
            Xfull(i, :) = Pt / cntLinks * Sel(i, :);
%             [~, Seli] = sort(Xi_ms2ms(i, :), 'descend');
%             Xfull(i, Seli(1:4)) = Pt / 4;
        end
        
        p_hat_full = pos_est_coop(Xfull, p_hat, J0_cell, Dist_ms2ms, Xi_ms2ms, Phi_ms2ms, ...
                  Pt_bs * Xi_ms2bs, Z_bs, Pos_bs);
        serr = zeros(1, Nms);
        for i_ms = 1:Nms
            serr(i_ms) = norm(p_hat_full(i_ms, :) - Pos_ms(i_ms, :))^2;
        end
		Mse_full(i_topo, i_sigma) = mean(serr); 
        Mse_ms1_full(i_topo, i_sigma) = serr(1);
        
        numLinks(i_sigma, 4) = length(find(Xfull > 0)) / 2;
		% --
        
        speb0_temp = 0;
        for k = 1:Nms
            speb0_temp = speb0_temp + trace(inv(J0_cell{k}));
        end
        Speb0(i_topo, i_sigma) = speb0_temp / Nms;
        
        X_full_naive = ones(Nms) / (Nms - 1) .* (ones(Nms) - eye(Nms));
        [Fspeb_full, ~] = obj_fun_np(X_full_naive, J0_cell, Xi_ms2ms, Phi_ms2ms, Vs(1));
        Speb1(i_topo, i_sigma) = mean(Fspeb_full);
        
        Speb_ms1_0(i_topo, i_sigma) = trace(inv(J0_cell{1}));
        Speb_ms1_1(i_topo, i_sigma) = Fspeb_full(1);

    end
    LinkNum(i_topo, :) = mean(numLinks, 1);
    
end

% % Cost (Objective values)
% Fc0 = repmat(Speb0, N_Vs, 1);
% Fc1 = repmat(Speb1, N_Vs, 1) + repmat(Vs(:), 1, N_sigma);
% Fc_g1 = Speb_g1 + repmat(Vs(:), 1, N_sigma) .* Pav_g1;
% Fc_g2 = Speb_g2 + repmat(Vs(:), 1, N_sigma) .* Pav_g2;
% Fc_g3 = Speb_g3 + repmat(Vs(:), 1, N_sigma) .* Pav_g3;
% Fc_ms1_g1 = Speb_ms1_g1 + repmat(Vs(:), 1, N_sigma) .* Pav_ms1_g1;
% Fc_ms1_g2 = Speb_ms1_g2 + repmat(Vs(:), 1, N_sigma) .* Pav_ms1_g2;
% Fc_ms1_g3 = Speb_ms1_g3 + repmat(Vs(:), 1, N_sigma) .* Pav_ms1_g3;
% Fc_ms1_0 = repmat(Speb_ms1_0, N_Vs, 1);
% Fc_ms1_1 = repmat(Speb_ms1_1, N_Vs, 1) + repmat(Vs(:), 1, N_sigma);
% 
% disp(Pav_g1);
% fprintf('No cooperation');
% disp(Speb0);
% fprintf('Full power cooperation');
% disp(Speb1);
% fprintf('Partial power cooperation');
% disp(Speb_g1);

%% Data processing
Spebs = [mean(Speb0)
         mean(Speb_g1)
         mean(Speb_g2)
         mean(Speb_g3)
         mean(Speb1)];
     
MSEs = [mean(Mse_nonc)
        mean(Mse_g1)
        mean(Mse_g2)
        mean(Mse_g3)
        mean(Mse_full)];

Speb_ue1 = [mean(Speb_ms1_0)
            mean(Speb_ms1_g1)
            mean(Speb_ms1_g2)
            mean(Speb_ms1_g3)
            mean(Speb_ms1_1)];
        
MSE_ue1 = [mean(Mse_ms1_nonc)
           mean(Mse_ms1_g1)
           mean(Mse_ms1_g2)
           mean(Mse_ms1_g3)
           mean(Mse_ms1_full)];
       
Pwrs = [mean(Pav_g1)
        mean(Pav_g2)
        mean(Pav_g3)
        ones(1, N_sigma);
        ];
Pwrs = Pwrs .* repmat(10 .^ (Sigma_dB / 10), 4, 1);
Pwrs_dB = 10 * log10(Pwrs);
%
my_markers = {'d', 's', '^', 'o', '<', '>', 'v', '+', 'p', '.', '*'}.';
my_colors = {'r', 'g', 'b', 'm', 'c', 'k'}.';

coop_scheme_order = [5 2 3 4];
    
Nsch = length(Scheme_str);
Scheme_speb = cell(Nsch, 1);
Scheme_mse = cell(Nsch, 1);
for i = 1:Nsch
    Scheme_speb{i} = sprintf('SPEB, %s', Scheme_str{i});
    Scheme_mse{i} = sprintf('MSE, %s', Scheme_str{i});
end

% % figure(1),
% if exist('fileHead', 'var') == 1
%     save_figure(hf1, fileHead, 'topology', {'fig','jpg'});
% end

% % Estimation error (SPEB + MSE) -------------------------------------------
% hf = figure(2);
% % SPEB --
% p_handle1 = plot(Sigma_dB, Spebs.', '--', 'LineWidth', 2, 'MarkerSize', 9);
% set(p_handle1, {'color'}, my_colors(1:size(Spebs, 1)));
% hold on
% % MLE --
% p_handle2 = plot(Sigma_dB, MSEs.', '-', 'LineWidth', 2, 'MarkerSize', 9);
% set(p_handle2, {'color'}, my_colors(1:size(MSEs, 1)));
% set(p_handle2, {'marker'}, my_markers(1:size(MSEs, 1)));
% hold off
% set(gca, 'FontSize',14);
% legend([Scheme_speb; Scheme_mse], 'FontSize', 12);
% xlabel('Power budget P_T (dB)');
% ylabel('MSE averaged over all users (m^2)');
% ylim([0 inf]);
% % grid on
% if exist('fileHead', 'var') == 1
%     save_figure(hf, fileHead, 'mse_all_ue', {'fig','jpg'});
% end
% 
% % Single User Estimation error (SPEB + MSE) -------------------------------
% hf = figure(3);
% % SPEB --
% p_handle1 = plot(Sigma_dB, Speb_ue1.', '--', 'LineWidth', 2, 'MarkerSize', 9);
% set(p_handle1, {'color'}, my_colors(1:size(Speb_ue1, 1)));
% hold on
% % MLE --
% p_handle2 = plot(Sigma_dB, MSE_ue1.', '-', 'LineWidth', 2, 'MarkerSize', 9);
% set(p_handle2, {'color'}, my_colors(1:size(MSE_ue1, 1)));
% set(p_handle2, {'marker'}, my_markers(1:size(MSE_ue1, 1)));
% hold off
% set(gca, 'FontSize',14);
% legend([Scheme_speb; Scheme_mse], 'FontSize', 12);
% xlabel('Power budget P_T (dB)');
% ylabel('MSE of User 1 (m^2)');
% ylim([0 inf]);
% % grid on
% if exist('fileHead', 'var') == 1
%     save_figure(hf, fileHead, 'mse_u1', {'fig','jpg'});
% end

% 
% % Cost value for all users ------------------------------------------------
% hf = figure(4);
% iVs = 1;
% % SPEB --
% plot(Sigma_dB, Fc0(iVs, :), 'k--', 'LineWidth', 1, 'MarkerSize', 8);
% hold on
% plot(Sigma_dB, Fc_g1(iVs, :), 'b--', 'LineWidth', 1, 'MarkerSize', 8);
% plot(Sigma_dB, Fc_g2(iVs, :), 'c--', 'LineWidth', 1, 'MarkerSize', 8);
% plot(Sigma_dB, Fc_g3(iVs, :), 'g--', 'LineWidth', 1, 'MarkerSize', 8);
% plot(Sigma_dB, Fc1(iVs, :), 'r--', 'LineWidth', 1, 'MarkerSize', 8);
% hold off
% set(p_handle, {'marker'}, my_markers(1:size(Mse_data, 1)));
% set(gca, 'FontSize',14);
% legend(Scheme_str, 'FontSize', 12);
% xlabel('Power budget P_T (dB)');
% ylabel('Objective values averaged over all users');
% ylim([0 inf]);
% % grid on
% if exist('fileHead', 'var') == 1
%     save_figure(hf, fileHead, 'cost_all_ue', {'fig','jpg'});
% end


% % Single User Cost -------------------------------------------------------
% hf = figure(5);
% iVs = 1;
% % SPEB --
% plot(Sigma_dB, Fc_ms1_0(iVs, :), 'k--', 'LineWidth', 1, 'MarkerSize', 8);
% hold on
% plot(Sigma_dB, Fc_ms1_g1(iVs, :), 'b--', 'LineWidth', 1, 'MarkerSize', 8);
% plot(Sigma_dB, Fc_ms1_g2(iVs, :), 'c--', 'LineWidth', 1, 'MarkerSize', 8);
% plot(Sigma_dB, Fc_ms1_g3(iVs, :), 'g--', 'LineWidth', 1, 'MarkerSize', 8);
% plot(Sigma_dB, Fc_ms1_1(iVs, :), 'r--', 'LineWidth', 1, 'MarkerSize', 8);
% hold off
% set(p_handle, {'marker'}, my_markers(1:size(Mse_data, 1)));
% set(gca, 'FontSize',14);
% legend(Scheme_str, 'FontSize', 12);
% xlabel('Power budget P_T (dB)');
% ylabel('Objective values of User 1');
% ylim([0 inf]);
% % grid on
% if exist('fileHead', 'var') == 1
%     save_figure(hf, fileHead, 'cost_u1', {'fig','jpg'});
% end

% % Power consumption of all users ------------------------------------------------
% hf = figure(6);
% iVs = 1;
% % Power --
% plot(Sigma_dB, 10 * log10(Pav_g1(iVs, :)), 'bd-', 'LineWidth', 1, 'MarkerSize', 8); hold on
% plot(Sigma_dB, 10 * log10(Pav_g2(iVs, :)), 'cs-', 'LineWidth', 1, 'MarkerSize', 8);
% plot(Sigma_dB, 10 * log10(Pav_g3(iVs, :)), 'g^-', 'LineWidth', 1, 'MarkerSize', 8);
% hold off
% set(gca, 'FontSize',14);
% legend({sprintf('Game1, V=%4.2f', Vs(iVs)), ...
%         sprintf('Game2, V=%4.2f', Vs(iVs)), ...
%         sprintf('Game3, V=%4.2f', Vs(iVs)), ...
%         }, 'FontSize', 12);
% xlabel('Power budget P_T (dB)');
% ylabel('Power consumed / power budget (dB)');
% ylim([-inf 0]);
% % grid on
% if exist('fileHead', 'var') == 1
%     save_figure(hf, fileHead, 'power_all_ue', {'fig','jpg'});
% end

% Estimation error MSE and power tradeoffs -------------------------------------------
hf = figure(7);
% % SPEB --
% p_handle1 = plot(Pwrs_dB, Spebs(2:end, :), '--', 'LineWidth', 2, 'MarkerSize', 9);
% set(p_handle1, {'color'}, my_colors(1:size(Spebs, 1) - 1));
% hold on
% MLE --
plot([0 Pwrs_dB(end, 2:end)].', MSEs(1, :), 'k--', 'LineWidth', 2, 'MarkerSize', 9);
hold on
p_handle2 = plot(Pwrs_dB.', MSEs(coop_scheme_order, :).', '-', 'LineWidth', 2, 'MarkerSize', 9);
set(p_handle2, {'color'}, my_colors(1:size(MSEs, 1) - 1));
set(p_handle2, {'marker'}, my_markers(1:size(MSEs, 1) - 1));
hold off
set(gca, 'FontSize',14);
legend([Scheme_mse(1); Scheme_mse(coop_scheme_order)], 'FontSize', 12, 'location', 'southwest');
xlabel('Actual power consumption (dBm)');
ylabel('MSE averaged over all agents (m^2)');
ylim([0 inf]);
% grid on
if exist('fileHead', 'var') == 1
    save_figure(hf, fileHead, 'mse_pwr_tof', {'fig','jpg'});
end

% % Estimation error MSE and power tradeoffs for agent 1 -------------------------------------------
% hf = figure(8);
% % % SPEB --
% % p_handle1 = plot(Pwrs_dB, Spebs(2:end, :), '--', 'LineWidth', 2, 'MarkerSize', 9);
% % set(p_handle1, {'color'}, my_colors(1:size(Spebs, 1) - 1));
% % hold on
% % MLE --
% p_handle2 = plot(Pwrs_dB.', MSE_ue1(2:end, :).', '-', 'LineWidth', 2, 'MarkerSize', 9);
% set(p_handle2, {'color'}, my_colors(1:size(MSE_ue1, 1) - 1));
% set(p_handle2, {'marker'}, my_markers(1:size(MSE_ue1, 1) - 1));
% hold off
% set(gca, 'FontSize',14);
% % legend([Scheme_speb; Scheme_mse], 'FontSize', 12);
% legend([Scheme_mse(2:end)], 'FontSize', 12);
% xlabel('Actual power consumption (dBm)');
% ylabel('MSE averaged over all users (m^2)');
% ylim([0 inf]);
% % grid on
% if exist('fileHead', 'var') == 1
%     save_figure(hf, fileHead, 'mse_pwr_tof_ue1', {'fig','jpg'});
% end


% %% ***************************** Send emails *****************************
% mail='huawei.hkust@gmail.com';
% password='testbed2577';
% 
% subject='Simulation job done! Data attached!';
% 
% emailto={'eejtchen@gmail.com'};
% 
% % Set up Gmail SMTP service.
% setpref('Internet','E_mail',mail);
% setpref('Internet','SMTP_Server','smtp.gmail.com');
% setpref('Internet','SMTP_Username',mail);
% setpref('Internet','SMTP_Password',password);
% 
% % Gmail server.
% props=java.lang.System.getProperties;
% props.setProperty('mail.smtp.auth','true');
% props.setProperty('mail.smtp.socketFactory.class','javax.net.ssl.SSLSocketFactory');
% props.setProperty('mail.smtp.socketFactory.port','465');
% 
% sendmail(emailto,subject,'See attached',{dataFileName});
