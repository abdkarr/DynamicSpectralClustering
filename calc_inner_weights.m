function weights = calc_inner_weights(A, C)
%CALC_INNER_WEIGHTS - Calculates number of intra-community edges of communities
%in a community structure of a undirected graph. If the graph is weighted, it
%calculates the sum of weights of intra-community edges.
%
%   Inputs:
%       A - nxn adjacency matrix of the graph.
%       C - n dimensional vector of community assignments.
%
%   Outputs:
%       weights - K dimensional column vector of intra-community edge number, 
%       where K is the number of communities in C. If the graph is weighted, 
%       it's sum of intra-community edge weights.
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: calc_association.m

%   Author: Abdullah Karaaslanli
%   Address: Michigan State University, ECE
%   email: karaasl1@msu.edu
%   Website: http://www.abdkarr.github.io
%   Date: 13-Nov-2020; Last revision: 13-Nov-2020
%
%   Copyright (c) 2020, Abdullah Karaaslanli
%   All rights reserved.

[comm_ids, n_comms] = get_comm_ids_number(C);

weights = zeros(n_comms, 1);
for r=1:n_comms
    comm_nodes = C == comm_ids(r); % node in the rth community
    weights(r) = sum(sum(A(comm_nodes, comm_nodes)))/2;
end
