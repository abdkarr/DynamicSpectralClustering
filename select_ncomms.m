function [g, ncomms] = select_ncomms(A, C, K, twom)
%SELECT_NCOMMS - Select the best community structure from a set of community 
%structures to determine number of communities. The best is determined as the
%one that has the maximum sum of modularity and asymptotic surprise.
%    
%   Inputs:
%       A - nxn adjacency matrix of the graph.
%       C - nxk dimensional matrix of community assignments that describe the 
%       set of community structures. k is the cardinality of the set.
%       twom - Twice the number of edges in the graph.
%
%   Outputs:
%       g - The best community structure.
%
%   Other m-files required: calc_modularity.m, calc_asurprise.me
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: 

%   Author: Abdullah Karaaslanli
%   Address: Michigan State University, ECE
%   email: karaasl1@msu.edu
%   Website: http://www.abdkarr.github.io
%   Date: 15-Apr-2021; Last revision: 15-Apr-2021
%
%   Copyright (c) 2021, Abdullah Karaaslanli
%   All rights reserved.

modularities = calc_modularity(A, C);
asurprises = calc_asurprise(A, C);
criterion = twom*modularities + asurprises;
[~, i] = max(criterion);

g = C(:, i);
ncomms = K(i);

end
