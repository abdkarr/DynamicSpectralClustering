function modularities = calc_modularity(A, C)
%CALC_MODULARITY - Calculates modularity value of each community structure in a 
%set of community structures of an undirected graph as defined in [1].
%
%   Inputs:
%       A - nxn adjacency matrix of the graph.
%       C - nxk dimensional matrix of community assignments that describe the 
%       set of community structures. k is the cardinality of the set.
%
%   Outputs:
%       modularities - k dimensional column vector, whose ith entry is modularity
%       value of the community structure C(:, k).
%
%   Other m-files required: calc_inner_weights.m, calc_association.m
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: calc_bipartite_modularity.m
%
%   References:
%       [1] Newman, M. E. J., and M. Girvan. “Finding and Evaluating Community 
%           Structure in Networks.” Physical Review E, vol. 69, no. 2, Feb. 
%           2004, p. 026113. DOI.org (Crossref), doi:10.1103/PhysRevE.69.026113.

%   Author: Abdullah Karaaslanli
%   Address: Michigan State University, ECE
%   email: karaasl1@msu.edu
%   Website: http://www.abdkarr.github.io
%   Date: 13-Nov-2020; Last revision: 29-Dec-2020
%
%   Copyright (c) 2020, Abdullah Karaaslanli
%   All rights reserved.

k = size(C, 2); % number of community structures provided

twom = sum(A(:));
modularities = zeros(k, 1);

for c=1:k
    modularities(c) = (2*sum(calc_inner_weights(A, C(:,c))) ...
        - sum(calc_association(A, C(:,c)).^2)/twom)/twom;  
end


