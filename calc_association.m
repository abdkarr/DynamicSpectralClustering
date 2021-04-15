function associations = calc_association(A, C)
%CALC_ASSOCIATION - Calculates associations of communities in a community 
%structure of an undirected graph. Association of a community is the sum of its
%nodes' degrees (or strengths in case the graph is weighted).
%    
%   Inputs:
%       A - nxn adjacency matrix of the graph.
%       C - n dimensional vector of community assignments.
%
%   Outputs:
%       associations - K dimensional column vector of association of communities,
%       where K is the number of communities in C.
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: calc_inner_weights.m

%   Author: Abdullah Karaaslanli
%   Address: Michigan State University, ECE
%   email: karaasl1@msu.edu
%   Website: http://www.abdkarr.github.io
%   Date: 13-Nov-2020; Last revision: 13-Nov-2020
%
%   Copyright (c) 2020, Abdullah Karaaslanli
%   All rights reserved.

[comm_ids, n_comms] = get_comm_ids_number(C);
d = sum(A, 2); % degree vector

associations = zeros(n_comms, 1);
for r=1:n_comms
    associations(r) = sum(d(C == comm_ids(r)));
end

end
