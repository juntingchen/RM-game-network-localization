function f = ls_coop_map(p, p0, J0, Lambda_c, Z_c, Pos_hat)
% f = ls_coop(p, Gamma, Z, Pos_bs, Lambda_c, Z_c, Pos_hat)

Nms = length(Z_c);

%f = 0;
% for i = 1:Nbs
%     f = f + Gamma(i) * (Z_bs(i) - norm(p - Pos_bs(i, :)))^2;
% end

f = (p - p0) * J0 * (p - p0).';

for j = 1:Nms
    f = f + Lambda_c(j) * (Z_c(j) - norm(p - Pos_hat(j, :)))^2;
end

% EXAMPLE on 3 users
% f = Gamma(1) * (Z_bs(1) - norm(p - Pos_bs(1, :)))^2 ...
%     + Gamma(2) * (Z_bs(2) - norm(p - Pos_bs(2, :)))^2 ...
%     + Gamma(3) * (Z_bs(3) - norm(p - Pos_bs(3, :)))^2 ...
%     + Lambda_c(1) * (Z_c(1) - norm(p - Pos_hat(1, :)))^2 ...
%     + Lambda_c(2) * (Z_c(2) - norm(p - Pos_hat(2, :)))^2;

    