function [comm_ids, n_comms] = get_comm_ids_number(c)
%GET_COMM_IDS_NUMBER - Get number of communities and ids of communities from
%an community assignment vector c.

%   Author: Abdullah Karaaslanli
%   Address: Michigan State University, ECE
%   email: karaasl1@msu.edu
%   Website: http://www.abdkarr.github.io
%   Date: 15-Apr-2021; Last revision: 15-Apr-2021
%
%   Copyright (c) 2021, Abdullah Karaaslanli
%   All rights reserved.

comm_ids = unique(c);
n_comms = length(comm_ids);

end
