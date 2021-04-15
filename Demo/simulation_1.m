%SIMULATION_1 - Detect communities of a simulated network generated using MLGM
%benchmark. The simulated network has transition probability(set to 0.9) and
%number of communities(set to 4) stationary across time. In the paper, number
%of nodes is set 128, mixing coefficient varies and minimum and maximum
%degrees of power law distribution respectively set to 8 and 16. Exponent of
%power-law distribution is set to -2.5. 
%
%Make sure MLGM (https://github.com/MultilayerGM/MultilayerGM-MATLAB) is
%downloaded and added to the path before running this script.
%
%   Other m-files required: none
%   MAT-files required: none
%
%   See also:

%   Author: Abdullah Karaaslanli
%   Address: Michigan State University, ECE
%   email: karaasl1@msu.edu
%   Website: http://www.abdkarr.github.io
%   Date: 30-Dec-2020; Last revision: 30-Dec-2020
%
%   Copyright (c) 2020, Abdullah Karaaslanli

n_runs = 100;

%% Set benchmark parameters
n_times = 15;
n_nodes = 128;

n_comms = 4;
% variation in community sizes, large mean eqauls sized communities
theta = 100;

% exponent for power-law distribution
exponent = -2.5;
min_degree = 8;
max_degree = 16;

mix_coeff = 0.5;

trans_probs = 0.9;
L = TemporalDependencyMatrix(n_times, trans_probs);

%% Algorithm parameters
K = 2:10; % candidate set of number of communities
max_iter = 20;

%% Outputs
nmi_on = zeros(n_runs, n_times);
nmi_off = zeros(n_runs, n_times);

% number of communities
nc_on = zeros(n_runs, n_times);
nc_off = zeros(n_runs, n_times);

%% Experiments
rng(46, 'twister');
reverse_str = '';
for r=1:n_runs
    % generate the network
    [A, gt] = DirichletDCSBMBenchmark(n_nodes, n_times, 'o', L, ...
        'UpdateSteps', 1, 'theta', theta, 'communities', n_comms, ...
        'exponent', exponent, 'kmin', min_degree, 'kmax', max_degree, ...
        'mu', mix_coeff, 'maxreject', 100);    
    A = cellfun(@full, A, 'UniformOutput', false);
    
    % find communities
    g_on = dsc_online(A, K, max_iter);
    g_off = dsc_offline(A, K, max_iter);
    
    % calculate nmi and number of communities: ignores disconnected nodes during
    % calculations
    for t=1:n_times
        conn_nodes = sum(A{t}) ~= 0; 
        
        nmi_on(r, t) = calc_nmi(g_on{t}(conn_nodes), gt(conn_nodes, t));
        nc_on(r, t) = length(unique(g_on{t}(conn_nodes)));
        
        nmi_off(r, t) = calc_nmi(g_off{t}(conn_nodes), gt(conn_nodes, t));
        nc_off(r, t)=length(unique(g_off{t}(conn_nodes)));
    end
    
    % Display the progress
    msg = sprintf('Experiment %d/%d is done.', r, n_runs); 
    fprintf([reverse_str, msg]);
    reverse_str = repmat(sprintf('\b'), 1, length(msg));
end

%% Plot the results
f = figure();

% Enlarge the figure to have 2 side by side subplots
f.Position(1) = f.Position(1) - f.Position(3)/2;
f.Position(3) = 2*f.Position(3);

subplot(1,2,1);
plot(1:n_times, mean(nmi_on, 1), 'DisplayName', '$\mathrm{DSC}_{on}$', ...
    'LineWidth', 1, 'Marker', 'o', 'MarkerFaceColor', 'auto');
hold on;
plot(1:n_times, mean(nmi_off, 1), 'DisplayName', '$\mathrm{DSC}_{off}$', ...
    'LineWidth', 1, 'Marker', 'o', 'MarkerFaceColor', 'auto');
l = legend;
l.Interpreter = 'latex';
l.Location = 'southeast';
title('NMI as a function of time', 'Interpreter', 'latex');
xlabel('Time', 'Interpreter', 'latex');
ylabel('NMI', 'Interpreter', 'latex');
grid;

subplot(1,2,2);
plot(1:n_times, mean(nc_on, 1), 'DisplayName', '$\mathrm{DSC}_{on}$', ...
    'LineWidth', 1, 'Marker', 'o', 'MarkerFaceColor', 'auto');
hold on;
plot(1:n_times, mean(nc_off, 1), 'DisplayName', '$\mathrm{DSC}_{off}$', ...
    'LineWidth', 1, 'Marker', 'o', 'MarkerFaceColor', 'auto');
l = legend;
l.Interpreter = 'latex';
l.Location = 'southeast';
title('Estimated number of communities', 'Interpreter', 'latex');
xlabel('Time', 'Interpreter', 'latex');
ylabel('Number of communities', 'Interpreter', 'latex');
grid;