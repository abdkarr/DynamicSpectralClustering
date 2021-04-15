function [sparsities] = calc_sparsities(A, C, twom)
%CALC_SPARSITIES - Calculates intra- and inter-sparsity of the graph using given
%community structure.
%
%   Inputs:
%       A - nxn adjacency matrix of the graph.
%       C - n dimensional vector of community assignments.
%       twom - Twice the number of edges in the graph.
%
%   Outputs:
%       sparsities - 2 dimensional column vector of intra- and inter-community
%       sparsities. sparsities(1) is the intra-community and sparsities(2) is 
%       the inter-community sparsity.
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none

%   Author: Abdullah Karaaslanli
%   Address: Michigan State University, ECE
%   email: karaasl1@msu.edu
%   Website: http://www.abdkarr.github.io
%   Date: 5-May-2020; Last revision: 4-Feb-2021
%
%   Copyright (c) 2020, Abdullah Karaaslanli
%   All rights reserved.

comm_ids = unique(C);
comm_ids(comm_ids==-1) = [];
num_comms = length(comm_ids);
num_nodes = length(C);
    
n_in = 0; % maximum number of intra-community edges
m_in = 0; % number of intra-community edges

for c=1:num_comms
    c_nodes = C == comm_ids(c);
    c_num_nodes = nnz(c_nodes);
    
    m_in = m_in + sum(sum(A(c_nodes, c_nodes)));
    n_in = n_in + c_num_nodes*(c_num_nodes-1);
end

n_out = num_nodes*(num_nodes-1) - n_in;
m_out = twom - m_in;

sparsities = zeros(2, 1);
sparsities(1) = m_in/(n_in+eps);
sparsities(2) = m_out/(n_out+eps);
   
end

