function thetas = estimate_thetas(A, C)
%ESTIMATE_THETAS - Estimate intra- and inter-community edge parameter of a
%degree-corrected planted partition model given binary undirected network and 
%community assignments. Degree correction is done with observed degrees. 
%
%   Inputs:
%       A - nxn adjacency matrix of the graph.
%       C - n dimensional vector of community assignments.
%
%   Outputs:
%       thetas - 2 dimensional vector of intra- and inter-community edge 
%       parameters. First entry is that of intra-community and second entry is
%       that of inter-community.
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: estimate_transitions.m

%   Author: Abdullah Karaaslanli
%   Address: Michigan State University, ECE
%   email: karaasl1@msu.edu
%   Website: http://www.abdkarr.github.io
%   Date: 30-Dec-2020; Last revision: 30-Dec-2020
%
%   Copyright (c) 2020, Abdullah Karaaslanli
%   All rights reserved.

assigned_nodes = C ~= -1;

degrees = sum(A, 2); 
m = sum(degrees(assigned_nodes))/2; % number of edges

m_in = sum(calc_inner_weights(A, C));
assocs = calc_association(A, C);

% intra-community edge parameter
thetas(1) = 2*m_in/(sum(assocs.^2)-sum(degrees.^2)+eps);
% inter-community edge parameter
thetas(2) = 2*(m-m_in)/(4*m^2-sum(assocs.^2)+eps); 

