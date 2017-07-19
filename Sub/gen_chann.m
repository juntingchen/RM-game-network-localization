function [Xi_ms2bs, Xi_ms2ms, Snr_ms2ms] = gen_chann(Dist_ms2bs, Dist_ms2ms)
% Generate the channel quality for the ranging (i.e., the equivalent
% ranging coefficient (ERC)).
% 
% This is built from Wenhan's code: topology_generation.m
% Reference in DaiSheWin14: "Network navigation algorithms with power
% control", WCNC, 2014.
%
% However, I did not understand why in [DaiSheWin14], uniform disbtituion
% on small scale fading is assumed. I modified that part.

Nms = size(Dist_ms2ms, 1);
Nbs = size(Dist_ms2bs, 2);

% Parameter Setting -------------------------------------------------------
Kfactor_dB = 4.7;
Kfactor = 10^(Kfactor_dB / 10);

los_ratio = 0.05;
rice_nu = sqrt(Kfactor / (1 + Kfactor));
rice_s = sqrt(1 / 2 / (1 + Kfactor));

BW = 20e6;
N0 = - 168;
Noise_dBm = - 168 + 10 * log10(BW);
Noise_Figure_dB = 5;

fc = 5.25e9;

% % LTE model delay line model
tau = [0 50 120 200 230 500 1600 2300 5000]; % delay profile
RP = [-1 -1 -1 0 0 0 -3 -5 -7];     % Power profile

% Winner model: B1 scenario: typical urban micro-cell
tau_bs = [0 30 55 60 105 115 250 460];
RP_bs = [0 -10.5 -14.8 -13.6 -13.9 -17.8 -19.6 -31.4];
tau_bs = tau;
% Winner model: A1 indoor office / residential 
tau_ms = [0 10 25 50 65 75 75 115 115 145 195 350];
RP_ms = [0 -15.1 -13.5 -15.1 -19.2 -23.5 -18.3 -23.3 -29.1 -14.2 -21.6 -23.4];
tau_ms = tau;
% Anchor parameters
f1 = BW / 1e9;
eff_f1 = f1/2/sqrt(3);
firstPath_bs = 10^(RP_bs(1)/10)/sum(10.^(RP_bs/10));

tau_len_bs = length(tau_bs);
A = zeros(tau_len_bs);
B = zeros(tau_len_bs);
C = zeros(tau_len_bs);
for i = 1 : tau_len_bs
    for j = 1 : tau_len_bs
        A(i,j) = sinc_new(tau_bs(i)-tau_bs(j),f1,2);
        B(i,j) = sinc_new(tau_bs(i)-tau_bs(j),f1,1);
        C(i,j) = sinc_new(tau_bs(i)-tau_bs(j),f1,0);
    end
end
J = [A -B; -B' C];
J_inv = inv(J);
chi_bs = 1/J_inv(1,1);

% Agent parameters
f2 = BW / 1e9;
eff_f2 = f2/2/sqrt(3);
firstPath_ms = 10^(RP_ms(1)/10)/sum(10.^(RP_ms/10));

tau_len_ms = length(tau_ms);
A = zeros(tau_len_bs);
B = zeros(tau_len_bs);
C = zeros(tau_len_bs);
for i = 1 : tau_len_ms
    for j = 1 : tau_len_ms
        A(i,j) = sinc_new(tau_ms(i)-tau_ms(j),f2,2);
        B(i,j) = sinc_new(tau_ms(i)-tau_ms(j),f2,1);
        C(i,j) = sinc_new(tau_ms(i)-tau_ms(j),f2,0);
    end
end
J = [A -B; -B' C];
J_inv = inv(J);
chi_ms = 1/J_inv(1,1);

% ERC for anchor measurements ---------------------------------------------
% Winner II, B1 LOS scenario
Xi_ms2bs = zeros(Nms, Nbs);
for i_bs = 1:Nbs
    for i_ms = 1:Nms
        pl = 41.0 + 22.7 * log10(Dist_ms2bs(i_ms, i_bs)) + 20 * log10(fc/1e9/5);
        sf = randn * 3.1;

        SNR = 10^((- pl - sf - Noise_dBm - Noise_Figure_dB) / 10) * firstPath_bs ...
              * ((1 - los_ratio) * random('rician', rice_nu, rice_s)^2 + los_ratio);
              % The los_ratio parameter is to avoid deep fade
              % In practice, deep fade can be removed by taking one more
              % measurement (no additional power)
        Xi_ms2bs(i_ms, i_bs) = chi_bs * 8 * pi^2 * (eff_f1*1e9)^2 / (3e8)^2 * SNR;
    end
end

% ERC for agent measurements ----------------------------------------------
% Winner II, A1 LOS scenario
Xi_ms2ms = zeros(Nms);
Snr_ms2ms = zeros(Nms);
for i = 1:Nms
    for j = i + 1:Nms
        pl = 46.8 + 18.7 * log10(Dist_ms2ms(i, j)) + 20 * log10(fc/1e9/5);
        sf = randn * 3.1;
        smallfade = random('rician', rice_nu, rice_s)^2;
        
        SNR = 10^((- pl - sf - Noise_dBm - Noise_Figure_dB) / 10) * firstPath_ms ...
              * smallfade;
        Xi_ms2ms(i, j) = chi_ms * 8 * pi^2 * (eff_f2*1e9)^2 / (3e8)^2 * SNR;
        Xi_ms2ms(j, i) = Xi_ms2ms(i, j);
        
        % SNR1 = 10^((- pl - sf - Noise_dBm - Noise_Figure_dB) / 10) * smallfade;
        Snr_ms2ms(i, j) = SNR;
        Snr_ms2ms(j, i) = SNR;
        
    end
end
