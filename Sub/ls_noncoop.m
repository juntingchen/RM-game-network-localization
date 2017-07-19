function f = ls_noncoop(p, Gamma, Z, Pos)
% f = ls_noncoop(p, Gamma, Z, Pos)

Nbs = length(Z);
f = 0;
for i = 1:Nbs
    f = f + Gamma(i) * (Z(i) - norm(p - Pos(i, :)))^2;
end


% EXAMPLE: 3 ANCHOR CASE
% f = Gamma(1) * (Z(1) - norm(p - Pos(1, :)))^2 ...
%     + Gamma(2) * (Z(2) - norm(p - Pos(2, :)))^2 ...
%     + Gamma(3) * (Z(3) - norm(p - Pos(3, :)))^2;