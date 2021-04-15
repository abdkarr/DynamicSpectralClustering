function Z = gen_indicator_mat(g)
%GEN_INDICATOR_MATRIX - Generates community indicator matrix from a given
%community assignment.
%
%   Syntax:
%       Z = gen_indicator_mat(g) generates nxK community indicator matrix from 
%           a given community assignment, where n is the number of nodes and K 
%           is the number of communities in the given assignment.
%        
%   Inputs:
%       g - n dimensional vector. Each of its entry is the id of the community 
%       the corresponding node belongs to. Entry with id -1 indicates node
%       doesn't belong to any communities.
%
%   Outputs:
%       Z - nxK indicator matrix. Each row has only one 1 and the rest is zero.
%       The column index of nonzero entry of the row is community assignment of
%       the node. If there are nodes not belong to any communities, their rows
%       are 0.
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none

%   Author: Abdullah Karaaslanli
%   Address: Michigan State University
%   email: karaasl1@msu.edu
%   Website: abdkarr.github.io
%   Date: 18-Feb-2020; Last revision: 24-Feb-2020
%
%   Copyright (c) 2020, Abdullah Karaaslanli
%   All rights reserved.

comm_ids = unique(g);
comm_ids(comm_ids == -1) = [];
num_comms = length(comm_ids);
num_nodes = length(g);

Z = zeros(num_nodes, num_comms);

for c = 1:num_comms
    Z(g == comm_ids(c), c) = 1; 
end

end
