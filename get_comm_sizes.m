function sizes = get_comm_sizes(C)
%GET_COMM_SIZES - Returns number of nodes in each community of a community
%structure.
%
%   Inputs:
%       C - n dimensional vector of community assignments.
%
%   Outputs:
%       sizes - K dimension column vector of community sizes, where K is the 
%       number of communties in C.
%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: 

%   Author: Abdullah Karaaslanli
%   Address: Michigan State University, ECE
%   email: karaasl1@msu.edu
%   Website: http://www.abdkarr.github.io
%   Date: 13-Nov-2020; Last revision: 13-Nov-2020
%
%   Copyright (c) 2020, Abdullah Karaaslanli
%   All rights reserved.

[~, ~, ic] = unique(C);
sizes = accumarray(ic, 1);


