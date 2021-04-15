function [transitions, num_remained, num_changed] = estimate_transitions(g_curr, ...
    g_prev)
%ESTIMATE_TRANSITIONS - Estimate intra- and inter-community transition
%statistics from previous time to current time point.
%    
%   Inputs:
%       g_curr - n dimensional vector of community assignment of nodes at current
%       time point. Its entries must be positive integers or -1. Positive
%       integers are assumed to be community ids and -1 means the node does not
%       belong to any community.
%       g_prev - n dimensional vector of community assignment of nodes at
%       previous time point. Its entries must be positive integers or -1. 
%       Positive integers are assumed to community ids and -1 means the node 
%       does not belong to any community.
%
%   Outputs:
%       transitions - 2 dimensional vector of intra- and inter-community
%       transition statistics. First entry is that of intra-community and second
%       is that of inter-community. 
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: estimate_thetas.m

%   Author: Abdullah Karaaslanli
%   Address: Michigan State University, ECE
%   email: karaasl1@msu.edu
%   Website: http://www.abdkarr.github.io
%   Date: 20-Apr-2020; Last revision: 20-Apr-2020
%
%   Copyright (c) 2020, Abdullah Karaaslanli
%   All rights reserved.

% pair of nodes that are in the same community
X_curr = bsxfun(@eq, g_curr, g_curr');
X_prev = bsxfun(@eq, g_prev, g_prev');

% ignore nodes that don't belong to any communities
outliers_curr = g_curr == -1;
outliers_prev = g_prev == -1;
outliers = outliers_prev | outliers_curr;

X_curr(outliers, :) = [];
X_curr(:, outliers) = [];
X_prev(outliers, :) = [];
X_prev(:, outliers) = [];

% Just consider upper triangle part, since it is undirected
X_curr = X_curr(triu(true(size(X_curr)), 1));
X_prev = X_prev(triu(true(size(X_prev)), 1));

% number of node pairs keep staying in the same community 
num_remained = nnz(X_curr & X_prev); 

% number of node pairs in the same community at previous time point
num_intra_pairs_prev = nnz(X_prev); 

% number of node pairs changed communities
num_changed = nnz(X_curr) - num_remained;

% calculate total pair of nodes at previous time
num_nodes_prev = length(g_prev)-nnz(outliers);
num_pairs_prev = num_nodes_prev*(num_nodes_prev-1)/2;

transitions = zeros(1, 2);
transitions(1) = num_remained/(num_intra_pairs_prev+eps);
transitions(2) = num_changed/ ...
                 (num_pairs_prev-num_intra_pairs_prev+eps);
            
end





