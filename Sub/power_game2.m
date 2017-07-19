function Xne = power_game2(J0_cell, Xi_mat, Phi_mat, V_vec, X0)
% [X1, Fspeb, Fcost] = power_game2(J0_cell, Xi_mat, Phi_mat, V_vec, X0)
%
% OUTPUT PARAMETERS
%   Xne     Nms * Nms matrix
%   Fspeb   1 * Nms row vector
%   Fcost   1 * Nms row vector, the objective values (cost)
% 
% VERSION 2: Fixed the channel gain issue for asynchronous networks
% 

Xi_mat = Xi_mat * 4; % Fix the channel gain issue by factor 4


NumberOfUpperLayerIterations = 3;

Nms = size(Xi_mat, 1);

% Form tree subnetworks
LinkSel0 = zeros(size(Xi_mat));
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
        if Xi_mat(k, j) <= Gamma(k, j)  % Necessary condition
            LinkSel0(k, j) = 0;
        else
            LinkSel0(k, j) = 1;
        end
    end
end

Bigxi = Xi_mat ./ (1e-11 + sqrt(Gamma) + sqrt(Gamma.')).^2;
LinkSel = zeros(Nms);
for k = 1:Nms
    [~, sel] = max(Bigxi(k, :));
    LinkSel(sel, k) = 1;
    LinkSel(k, sel) = 1;
end

Graph = LinkSel0 .* LinkSel;
NodeDegs = sum(Graph);

intNodes = find(NodeDegs > 1);      % Set of internal nodes
N_intNodes = length(intNodes);      
leNeighbor = cell(1, Nms);   % Set of leaf neighbor nodes of each internal node
intNeighbor = cell(1, Nms);  % Set of internal neighbor nodes of each internal node

node_mask = zeros(1, Nms);          % Indicator whether the node is visited
node_mask(intNodes) = 1;            % INternal nodes found

ConnectingInternalNode = 0;
for i = 1:N_intNodes
    neighbors = find(Graph(intNodes(i), :) > 0);    % Set of neighbor nodes of internal node i
    node_mask(neighbors) = 1;       % neighbors of internal node found
    leNeighbor{i} = []; % Set of leaf neighbor nodes of internal node i
    intNeighbor{i} = [];   % Set of internal neighbor nodes of internal node i
    for j = 1:length(neighbors)
        if NodeDegs(neighbors(j)) == 1
            leNeighbor{i} = [leNeighbor{i}, neighbors(j)];
        else
            ConnectingInternalNode = 1;
            intNeighbor{i} = [intNeighbor{i}, neighbors(j)];
        end
    end
end

% Only the two-user subnetworks left
cnt = N_intNodes;
while sum(node_mask) < Nms
    node_id = find(node_mask == 0, 1);
    intNodes = [intNodes, node_id];
    
    cnt = cnt + 1;
    leNeighbor{cnt} = find(Graph(node_id, :) > 0);
    node_mask([node_id, leNeighbor{cnt}]) = 1;
end
N_intNodes = cnt;

% Initialization
if nargin < 5
    X0 = zeros(Nms);
    for k = 1:Nms
        neighbors = find(Graph(k, :) > 0);
        X0(k, neighbors) = 1 / length(neighbors);
    end
end

% Fcost = zeros(1, Nms);
% Fspeb = zeros(1, Nms);

% Upper layer: Competitive game (Fixed number of iterations)
if ConnectingInternalNode
    Nloop = NumberOfUpperLayerIterations;
else 
    Nloop = 1;
end
for i = 1:Nloop
    X1 = X0;
    for i_node = 1:N_intNodes   % The i_node(th) internal node
        k = intNodes(i_node);
        leNb_k = leNeighbor{i_node};
        intNb_k = intNeighbor{i_node};
        
        ksubnet = [k, leNb_k, intNb_k];
        
        Xi_mat_k = Xi_mat(ksubnet, ksubnet);
        Phi_mat_k = Phi_mat(ksubnet, ksubnet);
        V_vec_k = V_vec(ksubnet);
        J0_cell_k = J0_cell(ksubnet);
        inode_mask = [zeros(1, length(leNb_k) + 1), ones(1, length(intNb_k))].';
        vec_p = X0(ksubnet, k);
        vec_p = vec_p .* inode_mask;
        vec_x0 = X0(k, ksubnet);
        
        [X_k, F_k] = power_ne_internal2(J0_cell_k, Xi_mat_k, Phi_mat_k, V_vec_k(1), ...
                                        vec_x0, vec_p, inode_mask);
        %
        X1(k, ksubnet) = X_k(1, :);
        X1(leNb_k, k) = X_k(2:length(leNb_k)+1, 1);
        % Fcost([k, leNb_k]) = F_k(1:1+length(leNb_k));
        
    end
    X0 = X1;
end

Xne = X1;

% for k = 1:Nms
%     Fspeb(k) = Fcost(k) - V_vec(k) * sum(Xne(k, :));
% end

