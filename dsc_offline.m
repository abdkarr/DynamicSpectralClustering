function g = dsc_offline(A, K, max_iter)
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

d = cellfun(@(X) sum(X, 2), A, 'UniformOutput', false); % degree vectors
twom = cellfun(@sum, d);
d_inv = cellfun(@(v) v.^-0.5, d, 'UniformOutput', false); % inverse degrees
conn_nodes = cellfun(@(v) v~=0, d, 'UniformOutput', false); % connected nodes

% data structs to keep output
g = cell(1, num_times);
theta = zeros(max_iter, 2, num_times);
pii = zeros(max_iter, 2, num_times-1);
num_comms = zeros(1, num_times);

%% FIND COMMUNITIES
for iter = 1:max_iter
    stop_nmi = true;
    stop_params = true;
    
    for t=1:num_times    
        % construct modified adjacency matrix
        A_mod = A{t}(conn_nodes{t}, conn_nodes{t});
    
        if iter==1
            % At first iteration modifiy adjacency only with previous time
            if t>1
                Z_prev = gen_indicator_mat(g{t-1});
                Z_prev = Z_prev(conn_nodes{t}, :);
                
                sparsities = calc_sparsities(A_mod, g{t-1}(conn_nodes{t}), ...
                    twom(t));
                S = (sparsities(1)-sparsities(2))*eye(num_comms(t-1)) + ...
                    sparsities(2);
                
                A_mod = A_mod + rm_diag((Z_prev*S*Z_prev'));
            end
        else
            % At remaining iterations, both next and previous time points are
            % used to modufy adjaceny matrix
            
            beta_t = log(theta(iter-1, 1, t)+eps)-log(theta(iter-1, 2, t)+eps);
            if t==1
                Z_next = gen_indicator_mat(g{t+1});
                Z_next = Z_next(conn_nodes{t}, :);
               
                sparsities = calc_sparsities(A{t+1}(conn_nodes{t+1}, ...
                    conn_nodes{t+1}), g{t}(conn_nodes{t+1}), twom(t+1));
                
                Pii_next = (sparsities(1)*pii(iter-1,1,t) ...
                    -sparsities(2)*pii(iter-1,2,t))*eye(num_comms(t+1));
                
                A_mod = beta_t*A_mod + rm_diag(Z_next*(Pii_next)*Z_next');
            elseif (t>1) && (t<num_times)
                Z_prev = gen_indicator_mat(g{t-1});
                Z_prev = Z_prev(conn_nodes{t}, :);
                Z_next = gen_indicator_mat(g{t+1});
                Z_next = Z_next(conn_nodes{t}, :);
                
                sparsities = calc_sparsities(A_mod, g{t-1}(conn_nodes{t}), ...
                   twom(t));
                S_prev = (sparsities(1)-sparsities(2))*eye(num_comms(t-1)) + ...
                    sparsities(2);
                Pii_prev = (pii(iter-1,1,t-1)-pii(iter-1,2,t-1))*eye(num_comms(t-1)) + ...
                    pii(iter-1,2, t-1);
                
                sparsities = calc_sparsities(A{t+1}(conn_nodes{t+1}, ...
                    conn_nodes{t+1}), g{t}(conn_nodes{t+1}), twom(t+1));
                Pii_next = (sparsities(1)*pii(iter-1,1,t)-sparsities(2)*pii(iter-1,2,t)) * ...
                    eye(num_comms(t+1));
                
                A_mod = beta_t*A_mod + rm_diag(Z_next*(Pii_next)*Z_next' + ...
                    Z_prev*(S_prev.*Pii_prev)*Z_prev');
            else
                Z_prev = gen_indicator_mat(g{t-1});
                Z_prev = Z_prev(conn_nodes{t}, :);
                
                sparsities = calc_sparsities(A_mod, g{t-1}(conn_nodes{t}), ...
                   twom(t));
                S_prev = (sparsities(1)-sparsities(2))*eye(num_comms(t-1)) + ...
                    sparsities(2);
                Pii_prev = (pii(iter-1,1,t-1)-pii(iter-1,2,t-1))*eye(num_comms(t-1)) + ...
                    pii(iter-1,2, t-1);
                
                A_mod = beta_t*A_mod+rm_diag(Z_prev*(S_prev.*Pii_prev)*Z_prev');
            end
        end
        
        % normalize adjacency matrix
        A_mod = diag(d_inv{t}(conn_nodes{t}))*A_mod*diag(d_inv{t}(conn_nodes{t}));
        
        % make sure the modified adjacency matrix is symmetric
        A_mod = (A_mod + A_mod')/2;
        
        % find communities
        if iter > 1
            g_old = g{t};
        end
        
        g{t} = -1*ones(num_nodes, 1);
        g_tmp = spectral_clustering(A_mod, K);
    
        % select number of communities
        [g{t}(conn_nodes{t}), num_comms(t)] = select_ncomms(...
            A{t}(conn_nodes{t}, conn_nodes{t}), g_tmp, K, twom(t));
        
        % check if the community structure is changed significantly compared
        % previous time point.
        if (iter > 1)
            if (calc_nmi(g_old, g{t}) < 0.9)
                stop_nmi = stop_nmi && false;
            else
                stop_nmi = stop_nmi && true;
            end
        end
        
        % estimate parameters
        if t==1
            theta(iter, :, t) = estimate_thetas(A{t}, g{t});
        else
            theta(iter, :, t) = estimate_thetas(A{t}, g{t});
            pii(iter, :, t-1) = estimate_transitions(g{t}, g{t-1});
        end
        
        % check the difference between new parameters and old parameters
        if iter>1 && t>1
            done = any((abs(theta(1:iter-1, 1, t)/theta(iter, 1, t) - 1) < 1e-2) & ...
                   (abs(theta(1:iter-1, 2, t)/theta(iter, 2, t) - 1) < 1e-2) & ...
                   (abs(pii(1:iter-1, 1, t-1)/pii(iter, 1, t-1) - 1) < 1e-2) & ...
                   (abs(pii(1:iter-1, 2, t-1)/pii(iter, 2, t-1) - 1) < 1e-2));
               
            if done
                stop_params = stop_params & true;
            else
                stop_params = stop_params & false;
            end
               
        end
    end
    
    % If community structure or parameters are not changing anymore for all time
    % points: algorithm is done.
    if iter>1
        if stop_nmi 
            break;
        end
        
        if stop_params
            break;
        end
    end
end

end
