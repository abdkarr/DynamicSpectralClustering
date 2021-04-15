function g = dsc_online(A, K, max_iter)
%DETECT_COMMS_ONLINE - Detect communities with online dynamic spectral
%clustering. See [1] for details. 
%
%   Inputs:
%       A - T dimensional cell array of adjacency matrices of dynamic network.
%       Make sure that number of nodes at each time is the same. If there are
%       emerging and disappearing nodes over time in your datasets add those to
%       all time points as disconnected nodes and the algorithm handles them by
%       not assigning them to any community.
%       K -  The candidate set of number of communities in the graph. Number of
%       communities at each time is selected from this set as the value that
%       maximizes the sum of modularity and asymptotic surprise.
%       max_iter - (Optional) Maximum number of iterations. Default is 20. 
%
%   Outputs:
%       g - T dimensional cell array of the assingment vector. g{t}[i] it the
%       id of community node i belongs to at time t. Note that, there is no
%       correspondance between community ids at consecutive time points. If this
%       is desired one can post-process communities with some algorithms
%       developped for this purpose. Finally, if g{t}[i] is equal to -1, that
%       means ith node at time t is a disconnected node and not assigned to any
%       community.
%
%   Other m-files required: gen_indicator_mat.m, calc_sparsities.m,
%   spectral_clustering.m, select_ncomms.m, estimate_thetas.m,
%   estimate_transitions.m, calc_nmi.m
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: dsc_offline.m
%
%   References:
%       [1] Karaaslanl?, Abdullah, and Selin Aviyente. "Community Detection in 
%           Dynamic Networks: Equivalence Between Stochastic Blockmodels and 
%           Evolutionary Spectral Clustering." IEEE Transactions on Signal and 
%           Information Processing over Networks 7 (2021): 130-143.

%   Author: Abdullah Karaaslanli
%   Address: Michigan State University, ECE
%   email: karaasl1@msu.edu
%   Website: http://www.abdkarr.github.io
%   Date: 17-Feb-2020; Last revision: 17-Feb-2020
%
%   Copyright (c) 2020, Abdullah Karaaslanli
%   All rights reserved.

if nargin < 3
    max_iter = 20;
end

rm_diag = @(X) X - diag(diag(X)); % removes the diagonal from the matrix

num_times = length(A);
num_nodes = size(A{1}, 1);

% data structs to keep output
g = cell(1, num_times); % found communities
theta = zeros(max_iter, 2, num_times); % estimated edge `probability values
pii = zeros(max_iter, 2, num_times-1); % estimated joint transition probabilities
num_comms = zeros(1, num_times);

%% FIND COMMUNITIES
for t=1:num_times
    
    % necessary parameters
    d = sum(A{t}, 2); % degree vector
    twom = sum(d); % two times number of edges
    conn_nodes = d~=0; % connected nodes
    
    %% INITIALIZING COMMUNITIES WITH PARAMETERS SET TO IDENTITY
    A_mod = A{t}(conn_nodes, conn_nodes);
    
    % modified adjacency with previous time community structure
    if t>1
        Z_prev = gen_indicator_mat(g{t-1});
        Z_prev = Z_prev(conn_nodes, :);
        
        sparsities = calc_sparsities(A_mod, g{t-1}(conn_nodes), twom);
        S = (sparsities(1)-sparsities(2))*eye(num_comms(t-1)) + sparsities(2);
        
        A_mod = A_mod + (rm_diag(Z_prev*S*Z_prev'));
    end
    
    % normalize the adjacency matrix
    A_mod = diag(d(conn_nodes).^-0.5)*A_mod*diag(d(conn_nodes).^-0.5);
    
    % make sure modified adjacency is symmetric
    A_mod = (A_mod + A_mod')/2;
    
    % find communities
    g{t} = -1*ones(num_nodes, 1);
    g_tmp = spectral_clustering(A_mod, K);
    
    % select number of communities
    [g{t}(conn_nodes), num_comms(t)] = select_ncomms(...
        A{t}(conn_nodes, conn_nodes), g_tmp, K, twom);
    
    %% COMMUNITIES INITIALIZED
    
    if t==1
        theta(1, :, t) = estimate_thetas(A{t}, g{t});
        continue;
    end
    
    %% ITERATIONS TO UPDATE COMMUNITIES WITH PARAMETER ESTIMATION
    for i = 1:max_iter
        
        % estimate parameters
        theta(i, :, t) = estimate_thetas(A{t}, g{t});
        pii(i, :, t-1) = estimate_transitions(g{t}, g{t-1});
        
        % if community structure at current and previous iterations are very 
        % similar or parameters didn't change much, break the iteration
        if i > 1
            done = calc_nmi(g_old(conn_nodes), g{t}(conn_nodes)) > 0.9;
            if done
                break;
            end
            
            done = any((abs(theta(1:i-1, 1, t)/theta(i, 1, t) - 1) < 1e-2) & ...
                   (abs(theta(1:i-1, 2, t)/theta(i, 2, t) - 1) < 1e-2) & ...
                   (abs(pii(1:i-1, 1, t-1)/pii(i, 1, t-1) - 1) < 1e-2) & ...
                   (abs(pii(1:i-1, 2, t-1)/pii(i, 2, t-1) - 1) < 1e-2));
            if done
                break;
            end
               
        end
        
        % construct modified adjacency
        beta_t = log(theta(i, 1, t) + eps) - log(theta(i, 2, t) + eps);
        Pi_t = (pii(i, 1, t-1)-pii(i, 2, t-1))*eye(size(Z_prev, 2)) ...
            + pii(i, 2, t-1);
        
        A_mod = beta_t*A{t}(conn_nodes, conn_nodes) + ...
            rm_diag(Z_prev*(S.*Pi_t)*Z_prev');
        
        % normalize the adjacency matrix
        A_mod = diag(d(conn_nodes).^-0.5)*A_mod*diag(d(conn_nodes).^-0.5);
        
        % make sure modified adjacency is symmetric
        A_mod = (A_mod + A_mod')/2;
        
        % find communities
        g_old = g{t};
        g{t} = -1*ones(num_nodes, 1);
        
        g_tmp = spectral_clustering(A_mod, K);
    
        % select number of communities
        [g{t}(conn_nodes), num_comms(t)] = select_ncomms(...
            A{t}(conn_nodes, conn_nodes), g_tmp, K, twom);
        
    end
    %% END OF ITERATIONS
end

end
