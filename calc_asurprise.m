function asurprises = calc_asurprise(A, C)
%CALC_ASURPRISE - Calculate asymptotic suprise of each community structure in a
%set of community structures of an undirected graph as defined in [1].
%
%   Inputs:
%       A - nxn adjacency matrix of the graph.
%       C - nxk dimensional matrix of community assignments that describe the 
%       set of community structures. k is the cardinality of the set.
%
%   Outputs:
%       asurprises - k dimensional vector, whose ith entry is asymptotic surprise
%       of the community structure C(:, k).
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: calc_modularity.m
%
%   References:
%       [1] Traag, Vincent A., Rodrigo Aldecoa, and J-C. Delvenne. "Detecting 
%           communities using asymptotical surprise." Physical Review E 92.2 
%           (2015): 022816.

%   Author: Abdullah Karaaslanli
%   Address: Michigan State University, ECE
%   email: karaasl1@msu.edu
%   Website: http://www.abdkarr.github.io
%   Date: 29-Dec-2020; Last revision: 29-Dec-2020
%
%   Copyright (c) 2020, Abdullah Karaaslanli
%   All rights reserved.

k = size(C, 2); % number of community structures provided

n = size(A, 1); % number of nodes
m = sum(A(:))/2; % number of edges
M = n*(n-1)/2; % number of node pairs

asurprises = zeros(k, 1);

for c=1:k
    sizes = get_comm_sizes(C(:, c));
    in_pairs = trace(sizes*(sizes-1)')/2; % number of intra-community node pairs
    m_in = sum(calc_inner_weights(A, C(:, c))); % number of intra-community edges
    
    q = m_in/m; % observed intra-community edge fraction
    q_avg = in_pairs/M; % expected intra-community edge fraction
    
    % KL divergence between q, q_avg
    asurprises(c) = m*(q*log((q+eps)/(q_avg+eps))+(1-q)*log((1-q+eps)/(1-q_avg+eps)));
end

