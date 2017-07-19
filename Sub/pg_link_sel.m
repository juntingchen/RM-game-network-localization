function [Graph, NodeDegs, GpDegs] = pg_link_sel(J0_cell, Xi_mat, Phi_mat, V_vec, L0)
% [Graph, NodeDegs, GpDegs] = pg_link_sel(J0_cell, Xi_mat, Phi_mat, V_vec)
%

Nms = size(Xi_mat, 1);

% Form tree subnetworks
% LinkSel0 = zeros(size(Xi_mat));
Gamma = zeros(Nms);
for k = 1:Nms
    Jk = J0_cell{k};
    Vk = V_vec(k);
    for j = 1:Nms
        if j == k
            continue;
        end
        t = Phi_mat(k, j);
        u = [cos(t), sin(t)].';
        Gamma(k, j) = Vk / (u' * Jk^(-2) * u);
% Comment: No need to perform necessary condition screening
%         if Xi_mat(k, j) <= Gamma(k, j)  % Necessary condition
%             LinkSel0(k, j) = 0;
%         else
%             LinkSel0(k, j) = 1;
%         end
    end
end

Bigxi = Xi_mat ./ (1e-11 + sqrt(Gamma) + sqrt(Gamma.')).^2;
LinkSel = zeros(Nms);
for k = 1:Nms
    [~, Isel] = sort(Bigxi(k, :), 'descend');
    sel = Isel(1: L0);
    LinkSel(sel, k) = 1;
    LinkSel(k, sel) = 1;
end

% Graph = LinkSel0 .* LinkSel;
Graph = LinkSel;

% Degrees of the nodes
NodeDegs = sum(Graph);

% Size of subnetworks
% Group detection via power of adjacent matrix
%   - Type A method: A is an adjacent matrix, [A]_{k,k} = 0, [A_{k,j} = 1
%   if k and j are connected. Property: [A^n]_{k,j} gives the number of
%   paths that nodes k and j are connected via exactly n edges.
%   
%   - Type B method: A is an adjacent matrix, [A]_{k,k} = 1, [A_{k,j} = 1
%   if k and j are connected. Property: [A^n]_{k,j} gives the number of
%   paths that nodes k and j are connected via AT MOST n edges.
A = eye(Nms) + LinkSel;

% To compute non-zero entries in Ap = A^(Nms - 1);
N = Nms - 1;
Ap = A;
for i = 1:floor(log2(N))
	Ap = Ap * Ap;
    Ap(Ap > 0) = 1;
end
for i = 1:N - round(2^floor(log2(N)))
    Ap = Ap * A;
    Ap(Ap > 0) = 1;
end

Gp = zeros(Nms);
node_record = zeros(1, Nms);
node_cnt = 0;
a_row_cnt = 1;
g_row_cnt = 0;
while node_cnt < Nms
    Arow = Ap(a_row_cnt, :);
    Arow_id = Arow > 0;
    node_record(Arow_id) = 1;
    if sum(node_record) > node_cnt
        % Found new group
        node_cnt = sum(node_record);
        g_row_cnt = g_row_cnt + 1;
        Gp(g_row_cnt, Arow_id) = 1;
    end
    a_row_cnt = a_row_cnt + 1;
end
Gp = Gp(1:g_row_cnt, :);
GpDegs = sum(Gp, 2);

% % -------------------------------------------------------------------------
% intNodes = find(NodeDegs > 1);      % Set of internal nodes
% N_intNodes = length(intNodes);      
% leNeighbor = cell(1, Nms);   % Set of leaf neighbor nodes of each internal node
% intNeighbor = cell(1, Nms);  % Set of internal neighbor nodes of each internal node
% 
% node_mask = zeros(1, Nms);          % Indicator whether the node is visited
% node_mask(intNodes) = 1;            % INternal nodes found
% 
% ConnectingInternalNode = 0;
% for i = 1:N_intNodes
%     neighbors = find(Graph(intNodes(i), :) > 0);    % Set of neighbor nodes of internal node i
%     node_mask(neighbors) = 1;       % neighbors of internal node found
%     leNeighbor{i} = []; % Set of leaf neighbor nodes of internal node i
%     intNeighbor{i} = [];   % Set of internal neighbor nodes of internal node i
%     for j = 1:length(neighbors)
%         if NodeDegs(neighbors(j)) == 1
%             leNeighbor{i} = [leNeighbor{i}, neighbors(j)];
%         else
%             ConnectingInternalNode = 1;
%             intNeighbor{i} = [intNeighbor{i}, neighbors(j)];
%         end
%     end
% end

